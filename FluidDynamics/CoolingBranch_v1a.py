import sys, os
sys.path.append('PressureDropAndHeatTransfer')
sys.path.append('../REFPROP')
import pandas as pd
import numpy as np
import matplotlib.pyplot as pl
from math import pi
from inspect import signature
from refprop import RefPropInterface
from time import time
from dPandHTC import *
from xml.dom import minidom
from FluidDynamicEquations import *

#use np arrays

import warnings
warnings.filterwarnings('ignore')

class CoolingBranch_v1a:
    refpropm = RefPropInterface().refpropm
    p=False
    def __init__(self, Fluid, Tsp, Tsc, Tsh, MF, xml, branch=None, varyValue=0, varyIndex=-1, eps=None, converge=None, param=None, dir=None):
        #config = {filename: , sheetname: }
        self.Fluid = Fluid
        self.setPointTemp = Tsp #C
        self.subCoolingTemp = Tsc #C
        self.allowedSuperHeatTemp = Tsh #C
        self.massFlow = MF #kg/s
#        self.heatFlow = HF #watts
#        self.HFf = HFf #?
        self.nargin = len(signature(CoolingBranch_v1a).parameters)
        self.branch=branch
        self.vary = varyIndex != -1
        if self.vary:
            self.dir = dir
            self.varyValue = float(varyValue)
            self.eps = float(eps)
            self.param = param
            self.converge = float(converge)
            self.dir = dir
            self.varyIndex = varyIndex
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
            "dL": 16
        }
        def getCol(N):
            parsed = []
            for j in range(len(params[N].getElementsByTagName("*"))):
                i = params[N].getElementsByTagName("*")[j]
                if i.tagName == "vary":
                    if self.vary:
                        # if self.dir == "up":
                        self.varyValue+=0.1
                    else:
                        self.varyIndex = i
                        self.varyValue = i.firstChild.data
                        self.eps = i.getAttribute("eps")
                        self.param = i.getAttribute("parameter")
                        self.dir = i.getAttribute("dir")
                        self.converge = i.getAttribute("converge")
                        self.vary = True

                    parsed.append(self.varyValue)
                    print("{}: {}".format(list(columns.keys())[j+1], self.varyValue))
                    # breakpoint()
                else:
                    parsed.append(i.firstChild.data)
            return np.array(list(map(float, parsed)))
        data = minidom.parse(xml)
        # breakpoint()
        params = data.getElementsByTagName("parameter") if not branch else data.getElementsByTagName("branch")[int(branch)].getElementsByTagName("parameter")
        # df = pd.read_excel(config['filename'], sheet_name=config['sheetname'])
        # getCol = lambda N: np.array(list(map(float, ["nan" if len(i.toxml()) == 8 else i.firstChild.data for i in params[N].getElementsByTagName("value")])))
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
        # print(self.ISOConductance)
        self.envTemperature=getCol(columns["envTemperature"]) #extract Tenv
        self.tubeWallThickness=getCol(columns["tubeWallThickness"])/1e3
        self.ISOThickness = getCol(columns["insulationThickness"])/1e3
        self.tubeThermalConductance = getCol(columns["tubeThermalConductance"])
        self.HTCAir = getCol(columns["HTCAir"])
        self.heatFlow=getCol(columns["heatFlow"])
        self.dL=getCol(columns["dL"])
    def start(self):
        self.initialize_arrays()
        self.fine_config()
        self.redefine()
    def initialize_arrays(self):
        tubeSurfaceArea = self.tubeSectionLength*np.pi*self.diameters
        self.heatFlux = self.heatFlow/tubeSurfaceArea #self.[float('nan')]*(len(self.heatSource))
        # for x in range(0, len(self.heatSource)):
        #     if self.heatSource[x] == 0:
        #         self.heatFlux[x] = 0
        #     elif all([np.isnan(i) for i in self.heatFlow]):
        #         self.heatFlux[x]=self.tubeWallThickness[x]/(self.tubeSectionLength[x]*np.pi*self.diameters[x]) #r_tube
        #     elif self.heatSource[x] == -1:
        #         self.heatFlux[x]=self.QPower[x]/(self.tubeSectionLength[x]*np.pi*self.diameters[x]) #Heatflux array
        #     else:
        #         self.heatFlux[x]=self.heatFlow[self.heatSource[x]-1]/(self.tubeSectionLength[x]*np.pi*self.diameters[x]) #Heatflux array
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

        # approximate dL to the closest value that allow for integer steps
        dLtm = self.tubeSectionLength/np.round(self.tubeSectionLength/self.dL)
        N = self.tubeSectionLength/dLtm


        # if sections will exchange heat, assume they have the same number of elements
        # if the section is only in thermal contact with air, HXNode == index of the section
        for i in range(len(N)):
            # print(self.HXNode[i])
            N[i] = N[self.HXNode[i]]
            dLtm[i] = self.tubeSectionLength[i]/N[i]

        # probably have to add the 1 as first element
        self.N = N
        self.SP = np.cumsum(N)
        self.SP = np.append([0],self.SP)
        self.SP = self.SP.astype(int)
        #self.SP+=1
        for i,n in enumerate(N):
            # print(n)
            # print(n, round(n), int(n))
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

        # is this necessary?
        # self.fineLength = self.fineLength[1:]
        self.fineHXNode = self.fineHXNode.astype(int)
        self.finePerimeterArea = np.pi*self.fineDiameter*self.finedLtm

        #
        # for x in range(len(self.tubeSectionLength)):
        #     if np.isnan(self.HXNode[x]):
        #         N = self.tubeSectionLength[x]/self.dL
        #     else:
        #         N = self.tubeSectionLength[self.HXNode[x]]/self.dL
        #     N=roundup(N)
        #     dLtm=self.tubeSectionLength[x]/N
        #     for y in range(len(self.fineLength), int(N)+len(self.fineLength)):
        #         self.fineLength.append(self.fineLength[y-1]+dLtm)
        #         self.fineDiameter.append(self.diameters[x])
        #         self.finePerimeterArea.append(np.pi*self.fineDiameter[-1]*dLtm)
        #         self.fineInclination.append(self.tubeOrientationAngle[x])
        #         self.fineRoughness.append(self.tubeRoughness[x])
        #         self.fineAppliedHeatFlux.append(self.heatFlux[x])
        #         self.fineHXNode.append(self.HXNode[x])
        #         self.fineHXFlowDir.append(self.HXFlowDir[x])
        #         self.fineHXConductance.append(self.HXConductance[x])
        #         self.fineInsulationConductance.append(self.ISOConductance[x])
        #         self.fineEnvTemperature.append(self.envTemperature[x])
        #     self.SP.append(len(self.fineLength))
        self.vaporQuality = np.zeros_like(self.fineLength) #[float('nan')]*len(self.fineLength)
        self.relativeMass = np.zeros_like(self.fineLength)
        self.dH = np.zeros_like(self.fineLength)
        self.wallTemperature = np.zeros_like(self.fineLength)
        self.fineHeatFlux = np.zeros_like(self.fineLength)
        self.State = ['']*len(self.fineLength)
        self.xia = np.zeros_like(self.fineLength)
        self.Gwavy = np.zeros_like(self.fineLength)
        self.Gwavy_xia = np.zeros_like(self.fineLength)
        self.Gstrat = np.zeros_like(self.fineLength)
        self.Gbub = np.zeros_like(self.fineLength)
        self.Gmist = np.zeros_like(self.fineLength)
        self.Gdry = np.zeros_like(self.fineLength)

        # self.rm = [float('nan')]*len(self.fineLength)
        # self.dH = [float('nan')]*len(self.fineLength)
        # self.wallTemperature = [float('nan')]*len(self.fineLength)
        # self.fineHeatFlux = [float('nan')]*len(self.fineLength)
        # self.State = ['']*len(self.fineLength)
    def redefine(self):
        self.fineHXNode_fCopy=self.fineHXNode #Redefine fine HX node arrangement taking flowdirection into account.
        # print(self.fineHXNode);
        # print(self.SP)
        # print(self.fineHXNode)
        for i in range(len(self.SP)-1):
            if self.fineHXFlowDir[self.SP[i]] == 0:
                # for j in range(self.SP[i], self.SP[i+1]):
                # print(self.SP)
                # print(self.SP[i], self.SP[i+1])
                self.fineHXNode_fCopy[self.SP[i]:self.SP[i+1]] = range(self.SP[i], self.SP[i+1])
                # try:
                #     self.fineHXNode_fCopy[self.SP[i]:self.SP[i+1]] = range(self.SP[i], self.SP[i+1])
                # except:
                #     breakpoint()
                # print(len(self.fineHXNode_fCopy[self.SP[i]:self.SP[i+1]]), len(range(self.SP[i], self.SP[i+1])))
            elif self.fineHXFlowDir[self.SP[i]] == -1:
                self.fineHXNode_fCopy[self.SP[i]:self.SP[i+1]] = range(self.SP[self.fineHXNode[self.SP[i]]+1]-1, self.SP[self.fineHXNode[self.SP[i]]]-1, -1)
            else:
                self.fineHXNode_fCopy[self.SP[i]:self.SP[i+1]] = range(self.SP[self.fineHXNode[self.SP[i]]], self.SP[self.fineHXNode[self.SP[i]]+1], 1)

        self.fineHXNode = self.fineHXNode_fCopy

        self.dP = np.ones_like(self.fineLength)
        self.HTC = np.ones_like(self.fineLength)
        _pressure = self.refpropm('P','T',self.setPointTemp+self.subCoolingTemp+273.15,'Q',0,self.Fluid)*1e-2;
        _enthalpy = self.refpropm('H','T',self.setPointTemp+self.subCoolingTemp+273.15,'Q',0,self.Fluid);
        self.P = np.ones_like(self.fineLength)*_pressure
        # subcooling temperature has to be set to 0 if CO2 is in saturation
        self.T = np.ones_like(self.fineLength)*(self.setPointTemp+self.subCoolingTemp)
        self.H = np.ones_like(self.fineLength)*_enthalpy

        # This is the set point information, which is fixed
        _setPointPressure = self.refpropm('P','T',self.setPointTemp+273.15,'Q',0,self.Fluid)*1e-2
        _setPointEnthalpy = self.refpropm('H','T',self.setPointTemp+273.15,'Q',0,self.Fluid)
        self.T[-1]=self.setPointTemp
        self.P[-1]=_setPointPressure
        self.H[-1]=_setPointEnthalpy

        self.fineMassFlux = self.massFlow/(0.25*np.pi*self.fineDiameter*self.fineDiameter)

        self.fineEnvHeatFlux = np.zeros_like(self.fineLength)
        self.fineHXHeatFlux = np.zeros_like(self.fineLength)

    def run(self, prt=True):
        format = lambda x: "'"+x[0]+"', "+",".join(map(str, x[1:]))
        itt=0;
        converge=10000;
        HXstart=0;
        convlimit=40;
        conv_repeat=0;
        conv_repeat_limit=3;
        self.HFLXf_hx = [float('nan')]*len(self.fineLength)
        self.HFf_hx = [0]*len(self.fineLength)
        self.Hconv = []
        self.HTCconv = []
        ittstop = 400
        avgHTC=0
        while (abs(converge)>convlimit or conv_repeat<conv_repeat_limit+1) and itt<ittstop:
            foo = True
            itt+=1
            if prt:
                print("\nIteration Round: {}".format(itt))
                print("Iteration offset: {} (Stops at {})".format(converge, convlimit))
                print("Enthalpy: {} J/kg".format(max(self.H)))
                print("Iteration offset: {} (Stops at {})".format(conv_repeat, conv_repeat_limit))
            Pprev = np.copy(self.P)
            Tprev = np.copy(self.T)
            HTCprev = np.copy(self.HTC)
            Hprev = np.copy(self.H)
            if converge < 1000:
                HXstart = 1
            times = {1: 0, 2: 0, 3: 0, 4: 0}
            startTime = time()
            columns, rows = os.get_terminal_size()
            for x in range(len(self.fineLength)-3,-1,-1): #2:length(Lf)
                #used to be range(len...-2, -1, -1)
                elapsed = time()-startTime
                total = elapsed*((len(self.fineLength)-2)-x)/(len(self.fineLength)-2)
                remaining = total-elapsed
                prefix = "\r{}/{} ({}%) [Elapsed time: {} Remaining time: {}] ".format(len(self.fineLength)-(x+2), len(self.fineLength)-2, round((len(self.fineLength)-x)/len(self.fineLength)*100), round(elapsed), round(remaining))
                progressbar_len = int(columns) - (len(prefix) + 2)
                progressbar = "".join(["#" if j <= (len(self.fineLength)-x)*progressbar_len//len(self.fineLength) else "." for j in range(progressbar_len)])
                #sys.stdout.write(prefix+"["+progressbar+"]");
                #sys.stdout.flush()
                # print(len(self.fineLength))
                # print("{}/{}".format(x, len(self.fineLength)))
                start = time()
                # print(self.fineHXNode[x+1])
                # print(x+1)
                # print(self.fineHXNode[x], x)
                if self.fineHXNode[x] == x:
                    # print(x)
                    #if in contact with self
                    avgHTC = (self.HTC[x]+self.HTC[x+1])/2
                    # print(self.HTC[x], self.HTC[x+1])
                    envResistance = 1/self.fineAirConductance[x] + self.fineInsulationThickness[x]/self.fineInsulationConductance[x]+self.fineTubeWallThickness[x]/self.fineTubeThermalConductance[x]+1/avgHTC

                    self.fineEnvHeatFlux[x] = (self.fineEnvTemperature[x]-self.T[x])/envResistance
                    # self.fineEnvHeatFlux[x] = 0
                    # if np.isnan(envResistance):
                    #     print(self.fineAirConductance[x], self.fineInsulationThickness[x], self.fineInsulationConductance[x], self.fineTubeWallThickness[x], self.fineTubeThermalConductance[x], avgHTC)
                    # print(avgHTC)
#_______________________________________
#                                       |
#              MODEL                    |
#                                       |
#              Air                      |
#              Insulation               |
#              Tube                     |
#              Fluid                    |
#                                       |
#_______________________________________|
                else:
                    nodeIndex = int(self.fineHXNode[x])
                    avgHTC = (self.HTC[x]+self.HTC[x+1])/2
                    avgNodeHTC = (self.HTC[nodeIndex]+self.HTC[nodeIndex+1])/2
                    hxResistance = 1/avgHTC+self.fineHXThickness[x]/self.fineHXConductance[x]+1/avgNodeHTC
                    self.fineHXHeatFlux[x] = (self.T[nodeIndex]-self.T[x])/hxResistance
                    # print(self.fineHXThickness[x], self.fineHXConductance[x], hxResistance, self)
                    # break;
                # print(self.fineEnvHeatFlux[x])
                # if np.isnan(self.fineHXHeatFlux[x]):
                #     # print(self.fineHXHeatFlux)
                #     break
                self.fineHeatFlux[x] = self.fineEnvHeatFlux[x]+self.fineHXHeatFlux[x]+self.fineAppliedHeatFlux[x]
                # print(self.fineEnvHeatFlux[x], self.fineHXHeatFlux[x], self.fineAppliedHeatFlux[x])
                if self.Fluid == 'CO2' and self.H[x] < 9e4: #avoid enthalpy to get into freezing area for cO2.
                    self.H[x] = 9e4

                newDP, newHTC, newVQ, newRM, newState, newT, newXia, newGwavy, newGwavy_xia, newGstrat, newGbub, newGmist, newGdry = dPandHTC(self.Fluid, round(self.P[x], 5), self.H[x], self.fineMassFlux[x], self.fineHeatFlux[x], self.fineDiameter[x], 0.25*pi*self.fineDiameter[x]**2, pi*self.fineDiameter[x], self.fineRoughness[x], self.fineInclination[x], self.allowedSuperHeatTemp, self.refpropm,x==316);
                # breakpoint()
                # print(newDP, newHTC, newVQ, newRM, newState, newT)
                # if x==316:
                #     print('debug1',newDP,newHTC,self.fineHeatFlux[x],self.P[x],round(self.P[x], 5),self.T[x],newT)
                # if(newDP > 100):
                #     print(newDP)
                #     print(format([self.Fluid, round(self.P[x], 5), self.H[x], self.fineMassFlux[x], self.fineHeatFlux[x], self.fineDiameter[x], 0.25*pi*self.fineDiameter[x]**2, pi*self.fineDiameter[x], self.fineRoughness[x], self.fineInclination[x], self.allowedSuperHeatTemp]))
                self.dP[x] = newDP
                # print(newDP)
                self.HTC[x] = newHTC
                # print(self.dP[x])
                # print("foo", self.HTC[x])
                self.vaporQuality[x] = newVQ
                self.relativeMass[x] = newRM
                self.State[x] = newState
                self.T[x] = newT
                # print(newT)
                self.xia[x] = float('nan') if not newXia else newXia.real
                self.Gwavy[x] = float('nan') if not newGwavy else newGwavy.real
                self.Gwavy_xia[x] = float('nan') if not newGwavy_xia else newGwavy_xia.real
                self.Gstrat[x] = float('nan') if not newGstrat else newGstrat.real
                self.Gbub[x] = float('nan') if not newGbub else newGbub.real
                self.Gmist[x] = float('nan') if not newGmist else newGmist.real
                self.Gdry[x] = float('nan') if not newGdry else newGdry.real
                # print(self.fineHeatFlux[x], self.fineDiameter[x], self.fineLength[x+1], self.fineLength[x])

                self.dH[x] = (self.fineHeatFlux[x]*pi*self.fineDiameter[x]*(self.fineLength[x]-self.fineLength[x+1]))/self.massFlow #calculate enthalpy difference, mass not mol
                # print(self.dH[x])
                if np.isnan(self.dH[x]):
                    print("ERR")
                    # print(self.Fluid, round(self.P[x], 5), self.H[x], self.fineMassFlux[x], self.fineHeatFlux[x], self.fineDiameter[x], 0.25*pi*self.fineDiameter[x]**2, pi*self.fineDiameter[x], self.fineRoughness[x], self.fineInclination[x], self.allowedSuperHeatTemp)
                    break
                # if newHTC == 0:
                #     print(self.Fluid, round(self.P[x+1], 5), self.H[x+1], self.fineMassFlux[x+1], self.fineHeatFlux[x+1], self.fineDiameter[x+1], 0.25*pi*self.fineDiameter[x+1]**2, pi*self.fineDiameter[x+1], self.fineRoughness[x+1], self.fineInclination[x+1], self.allowedSuperHeatTemp)
                #     # print(self.dH[x+1], self.fineEnvHeatFlux[x+1], self.fineHXHeatFlux[x+1], self.fineAppliedHeatFlux[x+1], self.T[x+1])
                #     print(self.HTC[x+1])
                #     break
                    # break
                avgHTC = (self.HTC[x]+self.HTC[x+1])/2
                partialResistance = self.fineTubeWallThickness[x+1]/self.fineTubeThermalConductance[x+1]+1/avgHTC
                self.wallTemperature[x] = self.T[x]+self.fineHeatFlux[x]*partialResistance #wall temperature
                self.H[x] = self.H[x+1]-self.dH[x] #calculate new enthalpy
                self.P[x] = self.P[x+1]+self.dP[x+1]*(self.fineLength[x+1]-self.fineLength[x]) #calculate new pressure
                # print(self.H[x+1])

            self.Hconv = np.array(self.H)-np.array(Hprev)
            self.HTCconv = np.array(self.HTC)-np.array(HTCprev)
            converge_prev = converge
            converge = max(self.Hconv) #use enthalpy to converge
            if abs(converge) < convlimit:
                conv_repeat+=1
            else:
                conv_repeat = 0
            # break
            #end while loop
    def plot(self):
        self.satTemperature = np.zeros_like(self.fineLength)
        for i in range(len(self.fineLength)):
            self.satTemperature[i] = self.refpropm('T','P',self.P[i]*1e2,'Q',self.vaporQuality[i], self.Fluid)-273.15

        #Intermittent to Annular Flow Transition Boundary
        fig, ax1 = pl.subplots(1)
        # fig, ax1 = pl.subplots()
        yax2 = ax1.twinx()
        ax1.plot(self.fineLength, self.T, 'g-', label='Temperature (Fluid)')
        ax1.plot(self.fineLength, self.wallTemperature, 'b-', label='Wall Temperature')
        ax1.plot(self.fineLength, self.satTemperature, 'c-', label='Saturated Temperature (Fluid)')
        yax2.plot(self.fineLength, self.P, 'r-', label='Pressure (bar)')
        ax1.legend()
        yax2.legend()
        ax1.set_xlabel('Length (m)')
        ax1.set_ylabel('Temperature (C)', color='g')
        yax2.set_ylabel('Pressure (bar)', color='b')
        fig2, ax2 = pl.subplots(1)

        ax2.axis(xmin=0, xmax=1)
        ax2.plot(self.vaporQuality, self.Gwavy, label='Gwavy')
        ax2.plot(self.vaporQuality, self.Gdry, label='Gdry')
        ax2.plot(self.vaporQuality[1:], self.fineMassFlux, '--', label='G')
        ax2.set_xlabel('Vapor Quality')
        ax2.set_ylabel('Mass Flux (kg/m^2s)')
        ax2.legend()
        # ax2.plot(self.fineLength, self.vaporQuality, 'r-')
        # ax2.set_xlabel('Length (m)')
        # ax2.set_ylabel('Vapor Quality (Fraction)')
        # fig3, ax3 = pl.subplots(1)
        # ax3.plot(self.vaporQuality, self.HTC, 'c-')
        # ax3.set_xlabel('Vapor Quality (Fraction)')
        # ax3.set_ylabel('Heat Transfer Coefficient kW/m^2K')
        # fig, apl.figure(1)
        # pl.
        # pl.scatter(self.P, self.T)
        # pl.xlabel("pressure")
        # pl.ylabel("temperature")
        # pl.figure(2)
        # pl.scatter(self.T, self.HTC)
        # pl.xlabel("temperature")
        # pl.ylabel("HTC")
        # pl.figure(3)
        # pl.scatter(x.HTC, x.fineHeatFlux)
        # pl.xlabel("HTC")
        # pl.ylabel("heat flux")
        pl.show(block=True)
        #P vs T (Tw)


if __name__ == "__main__":
    prefix = "../"
    filename = "CobraV1a_coupledring.xml" if len(sys.argv) <= 1 else sys.argv[1]
    path = prefix + filename
    MF = (0.7*16*20 + 0.7*4*9*2)/100
    x = CoolingBranch_v1a('CO2', -40, 0, 0, MF*1e-3, path)
    if x.vary:
        x.start()
        x.run(prt=False)
        cur = x.P[0]-x.P[-1]
        print("\tDelta Pressure: {}".format(cur))
        while cur > float(x.converge):
            # prevEps = x.eps
            # prevConverge = x.converge
            # prevParam = x.param
            # prevDir = x.dir
            x = CoolingBranch_v1a('CO2', -40, 0, 0, MF*1e-3, path, varyValue=x.varyValue, varyIndex=x.varyIndex, eps=x.eps, converge=x.converge, param=x.param, dir=x.dir)
            # x.eps = prevEps
            # x.converge = prevConverge
            # x.param = prevParam
            # x.dir = prevDir
            x.start()
            x.run(prt=False)
            cur = x.P[0]-x.P[-1]
            print("\tDelta Pressure: {}".format(cur))
    else:
        x.start()
        x.run()
    # print(x.P[0]-x.P[-1])
    # print("Completed. Plotting...")
    # print
    # x.plot()
