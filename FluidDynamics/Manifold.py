# from Node import *
# from SingleBranch import *
from xml.dom import minidom
import numpy as np
from refprop import RefPropInterface

class Manifold():
    
    refpropm = RefPropInterface().refpropm
    def __init__(self, base_path, root, initialMassFlow, initialVaporQuality=None, setPointTemp=None, subcoolingTemp=None, parent=None):

        from SingleBranch import SingleBranch
        self.SingleBranch = SingleBranch
        from Restrictor import Restrictor
        self.Restrictor = Restrictor
        
        self.Name = base_path+root
        self.Fluid = "CO2"
        self.setPointTemp = setPointTemp
        self.initialVaporQuality = initialVaporQuality or 0.01
        self.SubcoolingTemp = subcoolingTemp
        self.initialMassFlow = initialMassFlow
        self.branches = [] #array of manifolds (branches)
        self.xmls = []
        self.base_path = base_path
        self.root = root
        self.parent = parent
        self.series = None
        self.ST = None
        self.initialize()

    def initialize(self):
        self.branches = [(node.getAttribute('type'), node.getAttribute('loc')) for node in minidom.parse(self.base_path+self.root).getElementsByTagName("manifold")]
        self.series = minidom.parse(self.base_path+self.root).getElementsByTagName("components")[0].getAttribute("type") == "series";
        self.explore()

    def explore(self):
        for b in range(len(self.branches)):
            branch = self.branches[b]
            if branch[0] == 'manifold':
                self.branches[b] = Manifold(self.base_path, 
                                            branch[1], 
                                            self.initialMassFlow if self.series else self.initialMassFlow/len(self.branches), 
                                            initialVaporQuality=self.initialVaporQuality+(b+1)*(0.5-self.initialVaporQuality)/len(self.branches), 
                                            setPointTemp=self.setPointTemp, 
                                            subcoolingTemp=self.SubcoolingTemp,
                                            parent=self)
            else:
                self.branches[b] = self.SingleBranch(Fluid=self.Fluid, 
                                                     Tsp=self.setPointTemp, 
                                                     vq0 = self.initialVaporQuality if not self.series else self.initialVaporQuality+(b+1)*(0.5-self.initialVaporQuality)/len(self.branches), 
                                                     Tsc = self.SubcoolingTemp, 
                                                     MF = self.initialMassFlow if self.series else self.initialMassFlow/len(self.branches), 
                                                     xml = self.base_path+branch[1], 
                                                     parent=self)
    def updateMassFlow(self, newMF):
        self.initialMassFlow = newMF
        for branch in self.branches:
            if not self.series:
                branch.updateMassFlow(self.initialMassFlow/len(self.branches))
            else:
                branch.updateMassFlow(self.initialMassFlow)

    def minimize(self):
        N = len(self.branches)
        weights = np.full(N, self.initialMassFlow/N)
        for i in range(N):
            weights[i]=self.initialMassFlow*self.branches[i].getTotalAppliedHeat()/self.getTotalAppliedHeat()
            
        print("Total mass flow for {}: {} kg/s".format(" ".join(map(str, self.branches)), self.initialMassFlow))
        dPs = np.zeros_like(weights)
        epsilon = 0.005
        gamma = 0.0001
        counts = 0

        while dPs[0] == 0 or dPs.std() > epsilon:

            for b in range(len(self.branches)):
                branch = self.branches[b]
                branch.updateMassFlow(weights[b])
                branch.run(run=True)
                dPs[b] = branch.getDP()
            counts+=1
            print("Iteration {}:\n\tmass flows in each branch: {}\n\tdelta pressure in each branch: {}".format(counts, weights, dPs))
            avgDP = dPs.mean()
            diffs = dPs-np.full(N, avgDP)
            weights+=gamma*diffs
        return avgDP

    def getDP(self):
        return self.branches[-1].getFinalPressure() - self.branches[0].getStartPressure()
    def getFinalEnthalpy(self):
        return self.branches[-1].getFinalEnthalpy()
    def getFinalVaporQuality(self):
        return self.branches[-1].getFinalVaporQuality()
    def getInitialTemp(self):
        return self.branches[0].getInitialTemp()
    def getStartEnthalpy(self):
        return self.branches[0].getStartEnthalpy()
    def getStartPressure(self):
        return self.branches[0].getStartPressure()   
    def getFinalPressure(self):
        return self.branches[-1].getFinalPressure()  
    def getTotalAppliedHeat(self):
        retval = 0
        for b in self.branches:
            retval = retval + b.getTotalAppliedHeat()
        return retval
    
    def concat(self):

        # This series only work on two-phase.
        print('Solving concatenation')
        vaporQualities = np.linspace(self.initialVaporQuality, 0.5, len(self.branches)+1)
        for i in range(1,len(vaporQualities)):
            vaporQualities[i]=vaporQualities[i-1]+(0.5-self.initialVaporQuality)*self.branches[i].getTotalAppliedHeat()/self.getTotalAppliedHeat()

        temperatures = np.ones_like(vaporQualities)*self.setPointTemp
        print('Initial guesses')
        print('Vapor qualities: ', vaporQualities)
        print('Temperatures: ', temperatures)
        enthalpies = np.ones_like(vaporQualities)*self.refpropm('H','T',self.setPointTemp+273.15,'Q',self.initialVaporQuality,self.Fluid)
        print('Enthalpy: ', enthalpies)
        
        itt=0;
        converge=10000;
        convlimit=40;
        conv_repeat=0;
        conv_repeat_limit=2;
        ittstop = 400

        while (abs(converge)>convlimit or conv_repeat<conv_repeat_limit+1) and itt<ittstop:
            itt+=1
            print("\nManifold loop")
            print("Iteration Round: {}".format(itt))
            print("Iteration offset: {} (Stops at {})".format(converge, convlimit))
            print("Enthalpy: {} J/kg".format(max(enthalpies)))
            print("Iteration offset: {} (Stops at {})".format(conv_repeat, conv_repeat_limit))
            for ibranch in range(len(self.branches)-1, -1, -1):
                # these will remain unchanged
                self.branches[ibranch].setPointTemp = temperatures[ibranch+1]
                self.branches[ibranch].initialVaporQuality = vaporQualities[ibranch]
                # this will change
                self.branches[ibranch].finalVaporQualityGuess = vaporQualities[ibranch+1]
                # update
                print("Running Child")
                self.branches[ibranch].run(run = True)
                # update things that actually change
                vaporQualities[ibranch+1] = self.branches[ibranch].getFinalVaporQuality()            
                temperatures[ibranch] = self.branches[ibranch].getInitialTemp()

            print("Before updating")
            print('Vapor qualities: ', vaporQualities)
            print('Temperatures: ', temperatures)
            print('Enthalpy: ', enthalpies)
            prev_enthalpies = np.copy(enthalpies)
            for i,enthalpy in enumerate(enthalpies):
                enthalpies[i] = self.refpropm('H','T',temperatures[i]+273.15,'Q',vaporQualities[i],self.Fluid)
            print("Recalculating enthalpy")
            print('Vapor qualities: ', vaporQualities)
            print('Temperatures: ', temperatures)
            print('Enthalpy: ', enthalpies)

            total_dH = self.getStartEnthalpy() - self.refpropm('H','T',temperatures[0]+273.15,'Q',self.initialVaporQuality,self.Fluid);
            for i,enthalpy in enumerate(enthalpies):
                enthalpies[i] = enthalpy - total_dH
                vaporQualities[i] = self.refpropm('Q','T',temperatures[i]+273.15,'H',enthalpy,self.Fluid);
            print('Vapor qualities: ', vaporQualities)
            print('Temperatures: ', temperatures)
            print('Enthalpy: ', enthalpies)


            diff_enthalpies = np.abs(enthalpies-prev_enthalpies)
            converge = np.amax(diff_enthalpies) #use enthalpy to converge
            if abs(converge) < convlimit:
                conv_repeat+=1
            else:
                conv_repeat = 0

        return vaporQualities[-1]
       
    def __str__(self):
        return self.root
  
    def prettyprint(self, i=0, var=None):
        print(i*"\t|"+self.Name, end="    ")
        print(self.initialMassFlow if var == "IMF" else "")
        for b in self.branches:
            if b.__class__ == Manifold:
                b.prettyprint(i=i+1, var=var)
            else:
                print((i+1)*"\t|"+str(b), end="    ")
                print(b.massFlow if var == "IMF" else "")                
                
    def plot(self):
        print("Plotting manifold")
        for b in self.branches:
            b.plot()
            
    def run(self, ST=None, run=None):
#        print('RCLSA', self.series)
 
        ST = ST or self.setPointTemp
        
        if len(self.branches)==1:
            self.branches[0].run(run=True)
            print(self.getDP())
            return
        else:
            if self.series == True:
                print("Concatenating {}".format(" ".join(map(str, self.branches))))
                return self.concat()
            else:
                print("Minimizing {}".format(" ".join(map(str, self.branches))))
                return self.minimize()


