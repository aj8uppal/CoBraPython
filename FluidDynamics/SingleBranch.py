import sys, os
sys.path.append('../REFPROP')

import pandas as pd
import numpy as np
import matplotlib.pyplot as pl
import random
from math import pi
from inspect import signature
from time import time
from dPandHTC import *
from xml.dom import minidom
from FluidDynamicEquations import *
from Manifold import Manifold

#use np arrays
import warnings
warnings.filterwarnings('ignore')
np.set_printoptions(edgeitems=1000)

import runinfo

class SingleBranch(Manifold):
    
    p=False
    def __init__(self, Fluid, Tsp, vq0, vq0_point, MF, xml, parent=None, branch=None):
        
        # return
        self.xml = xml
        self.parent = parent
        self.Name = (xml.split('/')[-1])[:-4]
        self.Fluid = Fluid
        
        
        # This model the set point temperature of a PACL system
        self.setPointTemp = Tsp #C
        # These two model a warm nose
        self.initialVaporQuality = vq0 
        self.initialVaporQualitySector = vq0_point

        self.finalVaporQualityGuess = 0.5
        self.finalPressureGuess = 0
        self.finalUsePressureGuess = False
        
        self.branch=branch
        columns = {
            "section": 0,
            "hydraulicDiameter": 1, # in mm
            "tubeSectionLength": 2, # in m
            "tubeOrientationAngle": 3, # in degree
            "tubeRoughness": 4, # in um
            "HXNode": 5,
            "HXFlowDir": 6,
            "HXConductance": 7,
            "HXThickness": 8,
            "insulationConductance": 9,
            "envTemperature": 10,
            "tubeWallThickness": 11,
            "insulationThickness": 12,
            "tubeThermalConductance": 13,
            "HTCAir": 14,
            "heatFlow": 15,
            "dL": 16,
            "corrModel": 17 # 0: thome, 1: flex, 2: chisholm, 3: friedel, 4: thomex2 ?
        }
        def getCol(N):
            parsed = []
            for j in range(len(params[N].getElementsByTagName("*"))):
                i = params[N].getElementsByTagName("*")[j]
                parsed.append(i.firstChild.data)
            return np.array(list(map(float, parsed)))
        data = minidom.parse(xml)
        params = data.getElementsByTagName("parameter") if not branch else data.getElementsByTagName("branch")[int(branch)].getElementsByTagName("parameter")
        self.diameters=getCol(columns["hydraulicDiameter"])/1e3 # extract hydraulic diameter
        self.crossSectionArea=np.pi*self.diameters*self.diameters/4.
        self.tubeSectionLength=getCol(columns["tubeSectionLength"]) #extract lengths
        self.tubeOrientationAngle=getCol(columns["tubeOrientationAngle"]) #extract inclinations
        self.tubeRoughness=getCol(columns["tubeRoughness"])/1e6 #extract  roughnesses
        self.HXNode=np.array(list(map(int, getCol(columns["HXNode"])))) #extract heat exchange nodes
        self.HXFlowDir=getCol(columns["HXFlowDir"]) #extract heat exchanger flow direction
        self.HXConductance=getCol(columns["HXConductance"]) #extract heat exchanger conductance
        self.HXThickness=getCol(columns["HXThickness"])/1e3 #extract heat exchanger conductance
        self.ISOConductance=getCol(columns["insulationConductance"]) #extract insulation conductance
        self.envTemperature=getCol(columns["envTemperature"]) #extract Tenv
        self.tubeWallThickness=getCol(columns["tubeWallThickness"])/1e3
        self.ISOThickness = getCol(columns["insulationThickness"])/1e3
        self.tubeThermalConductance = getCol(columns["tubeThermalConductance"])
        self.HTCAir = getCol(columns["HTCAir"])
        self.heatFlow=getCol(columns["heatFlow"])
        self.dL=getCol(columns["dL"])
        self.corrModel=getCol(columns["corrModel"])
        self.massFlow = MF #kg/s
        if self.massFlow is None:
            self.massFlow = self.heatFlow.sum()*1e-5
    def updateMassFlow(self, newMF):
        self.massFlow = newMF
    def __eq__(self, other):
        return self.xml == other.xml
    def __repr__(self):
        return "<Single Branch {}>".format(self.xml)
    def __str__(self):
        return self.xml
    def run(self, ST=None, ivq=None, run=False):
        self.setPointTemp = ST if ST else self.setPointTemp
        self.initialVaporQuality = ivq if ivq else self.initialVaporQuality
        self.initialize_arrays()
        self.fine_config()
        self.redefine()
        if run:
            print('Child will run now')
            return self.main(prt=True)
    def initialize_arrays(self):
        tubeSurfaceArea = self.tubeSectionLength*np.pi*self.diameters
        self.heatFlux = self.heatFlow/tubeSurfaceArea #self.[float('nan')]*(len(self.heatSource))
    def fine_config(self):
        self.fineLength=[0]
        self.fineDiameter=[]
        self.fineInclination=[]
        self.fineHeatFlux=[]
        self.fineHXNode=[]
        self.fineHXFlowDir=[]
        self.fineHXConductance=[]
        self.fineAppliedHeatFlux = []
        self.fineInsulationConductance=[]
        self.fineEnvTemperature=[]
        self.finePerimeterArea = []
        self.fineRoughness = []
        self.fineAirConductance = []
        self.fineTubeWallThickness = []
        self.fineInsulationThickness = []
        self.fineTubeThermalConductance = []
        self.finedLtm = []
        self.fineHXThickness = []
        self.fineCorrModel = []

        # approximate dL to the closest value that allow for integer steps
        dLtm = self.tubeSectionLength/np.round(self.tubeSectionLength/self.dL)
        N = self.tubeSectionLength/dLtm

        # if sections will exchange heat, assume they have the same number of elements
        # if the section is only in thermal contact with air, HXNode == index of the section
        for i in range(len(N)):
            # print(self.HXNode[i])
            N[i] = N[self.HXNode[i]]
            dLtm[i] = self.tubeSectionLength[i]/N[i]

        self.N = N
        self.SP = np.cumsum(N)
        self.SP = np.append([0],self.SP)
        self.SP = self.SP.astype(int)
        
        self.initialVaporQualityFineSector = None
        if self.initialVaporQualitySector:
            self.initialVaporQualityFineSector = self.SP[int(self.initialVaporQualitySector)-1]

        for i,n in enumerate(N):
            n = int(round(n))
            self.fineLength = np.append(self.fineLength,
                np.linspace(self.fineLength[-1], self.fineLength[-1]+self.tubeSectionLength[i], n))
            self.fineDiameter = np.append(self.fineDiameter,[self.diameters[i]]*n)
            self.finedLtm = np.append(self.finedLtm, [dLtm[i]]*n)
            self.fineInclination = np.append(self.fineInclination, [self.tubeOrientationAngle[i]]*n)
            self.fineRoughness = np.append(self.fineRoughness, [self.tubeRoughness[i]]*n)
            self.fineAppliedHeatFlux = np.append(self.fineAppliedHeatFlux, [self.heatFlux[i]]*n)
            self.fineHXNode = np.append(self.fineHXNode, [int(self.HXNode[i])]*n)
            self.fineHXFlowDir = np.append(self.fineHXFlowDir, [self.HXFlowDir[i]]*n)
            self.fineHXConductance = np.append(self.fineHXConductance, [self.HXConductance[i]]*n)
            self.fineHXThickness = np.append(self.fineHXThickness, [self.HXThickness[i]]*n)
            self.fineInsulationConductance = np.append(self.fineInsulationConductance, [self.ISOConductance[i]]*n)
            self.fineEnvTemperature = np.append(self.fineEnvTemperature, [self.envTemperature[i]]*n)
            self.fineAirConductance = np.append(self.fineAirConductance,[self.HTCAir[i]]*n)
            self.fineTubeWallThickness = np.append(self.fineTubeWallThickness, [self.tubeWallThickness[i]]*n)
            self.fineInsulationThickness = np.append(self.fineInsulationThickness, [self.ISOThickness[i]]*n)
            self.fineTubeThermalConductance = np.append(self.fineTubeThermalConductance, [self.tubeThermalConductance[i]]*n)
            self.fineCorrModel = np.append(self.fineCorrModel, [self.corrModel[i]]*n)
 
        self.fineHeatFlux = np.zeros_like(self.fineAppliedHeatFlux)
        self.fineHXNode = self.fineHXNode.astype(int)
        self.finePerimeterArea = np.pi*self.fineDiameter*self.finedLtm
        self.vaporQuality = np.zeros_like(self.fineLength)
        self.relativeMass = np.zeros_like(self.fineLength)
        self.wallTemperature = np.zeros_like(self.fineLength)
        self.State = ['']*len(self.fineLength)
        self.xia = np.zeros_like(self.fineLength)
        self.Gwavy = np.zeros_like(self.fineLength)
        self.Gwavy_xia = np.zeros_like(self.fineLength)
        self.Gstrat = np.zeros_like(self.fineLength)
        self.Gbub = np.zeros_like(self.fineLength)
        self.Gmist = np.zeros_like(self.fineLength)
        self.Gdry = np.zeros_like(self.fineLength)


    def redefine(self):
        self.fineHXNode_fCopy=self.fineHXNode #Redefine fine HX node arrangement taking flowdirection into account.
        for i in range(len(self.SP)-1):
            if self.fineHXFlowDir[self.SP[i]] == 0:
                self.fineHXNode_fCopy[self.SP[i]:self.SP[i+1]] = range(self.SP[i], self.SP[i+1])
            elif self.fineHXFlowDir[self.SP[i]] == -1:
                self.fineHXNode_fCopy[self.SP[i]:self.SP[i+1]] = range(self.SP[self.fineHXNode[self.SP[i]]+1]-1, self.SP[self.fineHXNode[self.SP[i]]]-1, -1)
            else:
                self.fineHXNode_fCopy[self.SP[i]:self.SP[i+1]] = range(self.SP[self.fineHXNode[self.SP[i]]], self.SP[self.fineHXNode[self.SP[i]]+1], 1)

        self.fineHXNode = self.fineHXNode_fCopy
        self.fineMassFlux = self.massFlow/(0.25*np.pi*self.fineDiameter*self.fineDiameter)
        self.fineEnvHeatFlux = np.zeros_like(self.fineLength)
        self.fineHXHeatFlux = np.zeros_like(self.fineLength)

        self.HTC = 1000.*np.ones_like(self.fineLength)

        _pressure = self.refpropm('P','T',self.setPointTemp+273.15,'Q',self.initialVaporQuality,self.Fluid)*1e-2;
        _enthalpy = self.refpropm('H','T',self.setPointTemp+273.15,'Q',self.initialVaporQuality,self.Fluid);
        self.T = np.ones_like(self.fineLength)*self.setPointTemp
        self.P = np.ones_like(self.fineLength)*_pressure
        self.H = np.ones_like(self.fineLength)*_enthalpy

        # This is the set point information, which is fixed
        _setPointPressure = self.refpropm('P','T',self.setPointTemp+273.15,'Q',self.finalVaporQualityGuess,self.Fluid)*1e-2
        _setPointEnthalpy = self.refpropm('H','T',self.setPointTemp+273.15,'Q',self.finalVaporQualityGuess,self.Fluid)

        if self.finalUsePressureGuess:
            _setPointPressure = self.finalPressureGuess
            _setPointEnthalpy = self.refpropm('H','T',self.setPointTemp+273.15,'P',self.finalPressureGuess*1e2,self.Fluid)

        self.T[-1]=self.setPointTemp
        self.P[-1]=_setPointPressure
        self.H[-1]=_setPointEnthalpy

    def setFinalVaporQualityGuess(self, vq):
        self.finalVaporQualityGuess = vq
    def setFinalPressureGuess(self, pg):
        self.finalPressureGuess = pg

    def main(self, SH=None, ST=None, prt=False):
        # return 1
        format = lambda x: "'"+x[0]+"', "+",".join(map(str, x[1:]))
        itt=0;
        converge=10000;
        convlimit=40;
        conv_repeat=0;
        conv_repeat_limit=3;
        self.Hconv = []
        self.HTCconv = []
        ittstop = 400
        avgHTC=0
        while (abs(converge)>convlimit or conv_repeat<conv_repeat_limit+1) and itt<ittstop:
            itt+=1
            if prt:
                print("\n SB Iteration Round: {}".format(itt))
                print("SB Iteration offset: {} (Stops at {})".format(converge, convlimit))
                print("SB Enthalpy: {} J/kg".format(max(self.H)))
                print("SB Iteration offset: {} (Stops at {})".format(conv_repeat, conv_repeat_limit))
            Hprev = np.copy(self.H)

            for x in range(len(self.fineLength)-2,-1,-1):
                if self.fineHXNode[x] == x:
                    avgHTC = (self.HTC[x]+self.HTC[x+1])/2
                    envResistance = 1/self.fineAirConductance[x] + self.fineInsulationThickness[x]/self.fineInsulationConductance[x]+self.fineTubeWallThickness[x]/self.fineTubeThermalConductance[x]+1/avgHTC
                    self.fineEnvHeatFlux[x] = (self.fineEnvTemperature[x]-self.T[x+1])/envResistance
                else:
                    nodeIndex = int(self.fineHXNode[x])
                    avgHTC = (self.HTC[x]+self.HTC[x+1])/2
                    avgNodeHTC = (self.HTC[nodeIndex]+self.HTC[nodeIndex+1])/2
                    hxResistance = 1/avgHTC+self.fineHXThickness[x]/self.fineHXConductance[x]+1/avgNodeHTC
                    self.fineHXHeatFlux[x] = (self.T[nodeIndex]-self.T[x+1])/hxResistance
                self.fineHeatFlux[x] = self.fineEnvHeatFlux[x]+self.fineHXHeatFlux[x]+self.fineAppliedHeatFlux[x]
                newDP, newHTC, newVQ, newRM, newState, newT, newXia, newGwavy, newGwavy_xia, newGstrat, newGbub, newGmist, newGdry = dPandHTC(self.Fluid, self.P[x+1], self.H[x+1], self.fineMassFlux[x], self.fineHeatFlux[x], self.fineDiameter[x], 0.25*pi*self.fineDiameter[x]**2, pi*self.fineDiameter[x], self.fineRoughness[x], self.fineInclination[x], 0, self.refpropm, False, self.fineCorrModel[x]);
                try:
                    if np.isnan(newHTC):
                        raise ValueError("Oops, I will leave")
                except ValueError:
                    # if prt:
                    print('The code stepped away from its hypotheses')
                    print('-----------------------------------------')
                    print('Step ', x)
                    print('Length ', self.fineLength[x])
                    print('Vapor quality [old, new] ', self.vaporQuality[x+1], newVQ)
                    print('State [old, new] ', self.State[x+1], newState)
                    print('HTC [old, new] ', self.HTC[x+1], newHTC)
                    print('Temperature', self.T[x+1])
                    print('Pressure', self.P[x+1])
                    print('Enthalpy', self.H[x+1])
                    print('Environmental Heat Flux', self.fineEnvHeatFlux[x])
                    print('Heat Exchange Heat Flux', self.fineHXHeatFlux[x])
                    print('Applied Heat Flux', self.fineAppliedHeatFlux[x])
                    print('-----------------------------------------')
                    print('I will leave now. Bye!')
                    sys.exit(1)
                self.HTC[x+1] = newHTC
                self.vaporQuality[x+1] = newVQ
                self.relativeMass[x+1] = newRM
                self.State[x+1] = newState
                self.T[x+1] = newT
                self.Gwavy[x+1] = float('nan') if not newGwavy else newGwavy.real
                self.Gdry[x+1] = float('nan') if not newGdry else newGdry.real


                dH = (self.fineHeatFlux[x]*pi*self.fineDiameter[x]*(self.fineLength[x+1]-self.fineLength[x]))/self.massFlow #calculate enthalpy difference
                partialResistance = self.fineTubeWallThickness[x]/self.fineTubeThermalConductance[x]+1/self.HTC[x+1]
                self.wallTemperature[x+1] = self.T[x+1]+self.fineHeatFlux[x]*partialResistance #wall temperature
                self.H[x] = self.H[x+1]-dH #calculate new enthalpy
                self.P[x] = self.P[x+1]+newDP*(self.fineLength[x+1]-self.fineLength[x]) #calculate new pressure
                self.vaporQuality[x] = self.refpropm('Q','P',self.P[x]*1e2,'H',self.H[x],self.Fluid)
                self.T[x] = self.refpropm('T','P',self.P[x]*1e2,'H',self.H[x],self.Fluid)-273.15

            print('Calculating shift')
            total_dH = 0
            if self.initialVaporQualityFineSector:
                total_dH = self.H[self.initialVaporQualityFineSector] - self.refpropm('H','P',self.P[self.initialVaporQualityFineSector]*1e2,'Q',self.initialVaporQuality,self.Fluid)
            elif self.finalUsePressureGuess:
                total_dH = 0.0
            else:
                total_dH = self.H[0] - self.refpropm('H','P',self.P[0]*1e2,'Q',self.initialVaporQuality,self.Fluid)
            
            total_dH = random.uniform(0.8,1.0)*total_dH
            print("Shifting enthalpy: ", self.initialVaporQuality, self.H[0], total_dH)
            for i in range(len(self.H)):
                self.H[i] = self.H[i]-total_dH
                self.vaporQuality[i] = self.refpropm('Q','P',self.P[i]*1e2,'H',self.H[i],self.Fluid);

            self.Hconv = np.abs(self.H-Hprev)
            converge = np.amax(self.Hconv) #use enthalpy to converge

            if abs(converge) < convlimit:
                conv_repeat+=1
            else:
                conv_repeat = 0
                
        return self.getStartEnthalpy(), self.T[0]
    def getDP(self):
        return self.P[-1]-self.P[0]
    def getFinalEnthalpy(self):
        return self.H[-1]
    def getStartEnthalpy(self):
        return self.H[0]
    def getStartPressure(self):
        return self.P[0]
    def getFinalPressure(self):
        return self.P[-1]
    def getFinalVaporQuality(self):
        return self.vaporQuality[-1]
    def getInitialTemp(self):
        return self.T[0]
    def getTotalAppliedHeat(self):
        return self.heatFlow.sum()
    
    def plot(self):
        print("Plotting SB")
        print(self.T[self.initialVaporQualityFineSector])
        print(self.HTC[self.initialVaporQualityFineSector])
        self.satTemperature = np.zeros_like(self.fineLength)
        self.satPressure = np.zeros_like(self.fineLength)
        for i in range(len(self.fineLength)):
            self.satTemperature[i] = self.refpropm('T','P',self.P[i]*1e2,'H',self.H[i], self.Fluid)-273.15
            self.satPressure[i] = self.refpropm('P','T',self.T[i]+273.15,'H',self.H[i], self.Fluid)/1e2
            
        for i in range(len(self.fineLength)):
            if self.satPressure[i] < 0:
                self.satPressure[i] = self.satPressure[i-1]
                
        fig1, ax1 = pl.subplots(1)
        yax1 = ax1.twinx()
        ax1.plot(self.fineLength[1:], self.T[1:], 'g-', label='Temperature (Fluid)')
        #ax1.plot(self.fineLength[1:], self.satTemperature[1:], 'c-', label='Temperature (Saturation)')
        #ax1.plot(self.fineLength[1:-1], self.wallTemperature[1:-1], 'b-', label='Temperature (Wall)')
        yax1.plot(self.fineLength[1:], self.P[1:], 'r-', label='Pressure (Fluid)')
        yax1.plot(self.fineLength[1:], self.satPressure[1:], 'm-', label='Pressure (Saturation)')
        ax1.legend(loc='upper right')
        yax1.legend(loc='lower right')
        ax1.set_xlabel('Length (m)')
        ax1.set_ylabel('Temperature [C]', color='g')
        yax1.set_ylabel('Pressure [bar]', color='b')

        fig2, ax2 = pl.subplots(1)
        ax2.axis(xmin=0, xmax=self.vaporQuality[1:].max()*1.1, ymin=0, ymax=self.fineMassFlux.max()*1.1)
        ax2.plot(self.vaporQuality[1:], self.Gwavy[1:], label='Instability')
        ax2.fill_between(self.vaporQuality[1:], 0, self.Gwavy[1:], facecolor='blue', alpha=0.3)
        ax2.plot(self.vaporQuality[1:], self.Gdry[1:], label='Dry-out')
        ax2.fill_between(self.vaporQuality[1:], self.Gdry[1:], self.fineMassFlux.max()*2, facecolor='orange', alpha=0.3)
        ax2.plot(self.vaporQuality[1:], self.fineMassFlux, '--', label='G')
        ax2.set_xlabel('Vapor Quality')
        ax2.set_ylabel('Mass Flux (kg/m^2s)')
        ax2.legend()

        fig3, ax3 = pl.subplots(1)
        yax3 = ax3.twinx()
        ax3.plot(self.fineLength[1:], self.HTC[1:], 'g-', label='Heat Transfer Coef')
        yax3.plot(self.fineLength[1:], self.vaporQuality[1:], 'r-', label='Vapor Quality')
        ax3.legend(loc='upper left')
        yax3.legend(loc='lower right')
        ax3.set_xlabel('Length (m)')
        ax3.set_ylabel('HTC (W/m^2K)', color='g')
        yax3.set_ylabel('Vapor Quality', color='r')

        fig1.savefig('output/' + runinfo.runname + '_' + self.Name + '_PT.pdf')
        fig2.savefig('output/' + runinfo.runname + '_' + self.Name + '_Flow.pdf')
        fig3.savefig('output/' + runinfo.runname + '_' + self.Name + '_HTC.pdf')


        print('SB writing data')
        df_columns = ['Pressure', 'Temperature', 'Enthalpy', 'HTC', 'satPressure', 'satTemperature', 'wallTemperature', 'vaporQuality', 'heatFlux', 'Length', 'State']
        thermalData = []
        for i in range(len(self.P)-1):
            thermalData.append([self.P[i], self.T[i], self.H[i], self.HTC[i], self.satTemperature[i], self.satPressure[i], self.wallTemperature[i], self.vaporQuality[i], self.fineHeatFlux[i], self.fineLength[i], self.State[i]])
        thermal_df = pd.DataFrame(thermalData, columns=df_columns)
        thermal_df.to_hdf('output/' + runinfo.runname + '_' + self.Name + '.h5', key='fluid', mode='w')
        thermal_df.to_csv('output/' + runinfo.runname + '_' + self.Name + '.csv', mode='w')


        
