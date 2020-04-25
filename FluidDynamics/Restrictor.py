from xml.dom import minidom
import numpy as np
from refprop import RefPropInterface
from Manifold import Manifold

class Restrictor(Manifold):
    
    def __init__(self, Fluid, Tsp, vq0, Tsc, xml, parent=None, branch=None, dir=None):

        self.Name = (xml.split('/')[-1])[:-4]
        self.Fluid = Fluid
        self.dP = None
        self.setPointTemp = Tsp #C
        self.initialVaporQuality = vq0
        self.SubcoolingTemp = Tsc #C
        self.finalVaporQualityGuess = 0.5
        self.xml = xml
        self.parent = parent
        self.initialize()

    def initialize(self):
        data = minidom.parse(self.xml)
        params = data.getElementsByTagName("parameter")
        self.dP = params[0].getAttribute("PressureDrop")
        self.initialTemperature = self.setPointTemp
        self.finalTemperature = self.setPointTemp
        self.initialEnthalpy = 0
        self.finalEnthalpy = 0
        self.initialPressure = 0
        self.finalPressure = 0

    def getDP(self):
        return self.dP
    def getFinalEnthalpy(self):
        return self.finalEnthalpy
    def getFinalVaporQuality(self):
        return self.initialVaporQuality
    def getInitialTemp(self):
        return self.initialTemperature
    def getStartEnthalpy(self):
        return self.initialEnthalpy
    def getStartPressure(self):
        return self.initialPressure
    def getFinalPressure(self):
        return self.finalPressure 
    def getTotalAppliedHeat(self):
        return 0
                
    def plot(self):
        return
            
    def run(self, ST=None, run=None):
        
        if run is None or run is False:
            return 0

        
        self.finalTemperature = ST if ST else self.setPointTemp        
        
        finalPressure = 0
        if self.SubcoolingTemp:
                        
            self.initialTemperature = self.SubcoolingTemp
            convLimit = 40
            eps = 0.01
            self.intialPressure = self.refpropm('P','T',self.finalTemperature+273.15,'Q',self.finalVaporQualityGuess,self.Fluid)-self.dP
            
            while self.dH > convLimit:
                self.finalEnthalpy = self.refpropm('H','T',self.finalTemperature+273.15,'P',self.initialPressure+self.dP,self.Fluid)
                self.initialEnthalpy = self.refpropm('H','T',self.intialTemperature+273.15,'P',self.initialPressure,self.Fluid)
                self.dH = self.finalEnthalpy - self.initialEnthalpy
                self.initialPressure = self.initialPressure - eps
            self.finalPressure = self.initialPressure + self.dP
            
        else:
            self.finalPressure = self.refpropm('P','T',self.finalTemperature+273.15,'Q',self.initialVaporQuality,self.Fluid)
            self.finalEnthalpy = self.refpropm('H','T',self.finalTemperature+273.15,'Q',self.initialVaporQuality,self.Fluid)          
            self.initialEnthalpy = self.finalEnthalpy
            self.initialPressure = self.finalPressure - self.dP             
            self.initialTemperature = self.refpropm('T','P',self.initialPressure*1e2,'Q',self.initialVaporQuality,self.Fluid)
            
        return 0