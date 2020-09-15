import sys
sys.path.append('../REFPROP')

from xml.dom import minidom
import numpy as np
from refprop import RefPropInterface

#from multiprocessing import Pool
import multiprocessing
import multiprocessing.pool

class NoDaemonProcess(multiprocessing.Process):
    @property
    def daemon(self):
        return False

    @daemon.setter
    def daemon(self, value):
        pass

class NoDaemonContext(type(multiprocessing.get_context())):
    Process = NoDaemonProcess

class MyPool(multiprocessing.pool.Pool):
    def __init__(self, *args, **kwargs):
        kwargs['context'] = NoDaemonContext()
        super(MyPool, self).__init__(*args, **kwargs)
        
class Manifold():
    
    refpropm = RefPropInterface().refpropm
    def __init__(self, base_path, root, initialMassFlow, initialVaporQuality=None, initialVaporQualitySector = None, setPointTemp=None, parent=None):

        from SingleBranch import SingleBranch
        self.SingleBranch = SingleBranch
        
        self.Name = base_path+root
        self.Fluid = "CO2" # this has to be determined by the XML... update later.
        self.setPointTemp = setPointTemp
        self.initialVaporQuality = initialVaporQuality #put 1% just to be safe
        self.initialVaporQualitySector = initialVaporQualitySector
        self.initialMassFlow = initialMassFlow
        self.branches = [] #array of manifolds (branches)
        self.xmls = []
        self.base_path = base_path
        self.root = root
        self.parent = parent
        self.series = None
        self.initialize()

    def initialize(self):
        self.branches = [[node.getAttribute('type'), node.getAttribute('loc'), node.getAttribute('vq0'), node.getAttribute('vq0_point')] for node in minidom.parse(self.base_path+self.root).getElementsByTagName("manifold")]
        for branch in self.branches:
            branch[2] = float(branch[2])
            branch[3] = int(branch[3])
        self.series = minidom.parse(self.base_path+self.root).getElementsByTagName("components")[0].getAttribute("type") == "series";
        self.explore()
        
    def explore(self):
        for b in range(len(self.branches)):
            branch = self.branches[b]
            if branch[0] == 'manifold':
                self.branches[b] = Manifold(self.base_path, 
                                            branch[1], 
                                            self.initialMassFlow if self.series else self.initialMassFlow/len(self.branches), 
                                            0.0, # Manifold never have initial vapor quality, this is defined at the singlebranch level
                                            None,
                                            setPointTemp=self.setPointTemp, 
                                            parent=self)
            else:                
                self.branches[b] = self.SingleBranch(Fluid=self.Fluid, 
                                                     Tsp=self.setPointTemp, 
                                                     vq0 = branch[2],
                                                     vq0_point = branch[3] if branch[3]>=0 else None,                                                     
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

    def worker(self, branch):
        branch.run(run=True)
        return branch
    
    def minimize(self):
        N = len(self.branches)
        weights = np.full(N, self.initialMassFlow/N)
        for i in range(N):
            weights[i]=self.initialMassFlow*self.branches[i].getTotalAppliedHeat()/self.getTotalAppliedHeat()
            
        print("Total mass flow for {}: {} kg/s".format(" ".join(map(str, self.branches)), self.initialMassFlow))
        dPs = np.zeros_like(weights)
        epsilon = 0.1
        gamma = 0.0001
        counts = 0

        while dPs[0] == 0 or dPs.std() > epsilon:
            
            for b in range(len(self.branches)):
                branch = self.branches[b]
                branch.updateMassFlow(weights[b])
                # branch.run(run=True)

            p = MyPool(len(self.branches))
            self.branches = p.map(self.worker, self.branches)

            for b in range(len(self.branches)):
                branch = self.branches[b]
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

        print('Solving concatenation')

        # logic:
        #   solve each part of the series
        #   from nBranch-1
        #      if it is evaporating
        #            set setPointTemp to initial temp of the next one
        #            set the final guess to the initial vapor quality of the
        #      if not evaporating
        #            set setPointTemp to same temperature as the final branch
        #            set final guess to 0.
        
        
        for ibranch, branch in enumerate(self.branches):
            if branch.initialVaporQuality > 0:
                firstEvapBranch = ibranch
                firstVaporQuality = branch.initialVaporQuality
                break
        nEvapBranches = len(self.branches) - firstEvapBranch;
        
        vaporQualities = np.linspace(self.initialVaporQuality, 0.5, nEvapBranches+1)
        for i in range(firstEvapBranch+1,len(self.branches)):
            vaporQualities[i-firstEvapBranch]=vaporQualities[i-1-firstEvapBranch]+(0.5-self.branches[firstEvapBranch].initialVaporQuality)*self.branches[i-1].getTotalAppliedHeat()/self.getTotalAppliedHeat()
        vaporQualities = np.concatenate(([0]*firstEvapBranch,vaporQualities))
        temperatures = np.ones_like(vaporQualities)
        pressures = np.ones_like(vaporQualities)
        for i in range(len(self.branches)):
            temperatures[i+1] = self.branches[i].setPointTemp
        temperatures[0] = self.branches[0].setPointTemp #getInitialTemp() 
        for i in range(len(self.branches)+1):
            pressures[i] = self.refpropm('P','T',temperatures[i]+273.15,'Q',vaporQualities[i],self.Fluid)*1e-2;
        print('Initial guesses')
        print('Vapor qualities: ', vaporQualities)
        print('Temperatures: ', temperatures)

        itt=0;
        converge=10000;
        convlimit=100;
        conv_repeat=0;
        conv_repeat_limit=3;
        ittstop = 400

        while (abs(converge)>convlimit or conv_repeat<conv_repeat_limit+1) and itt<ittstop:
            itt+=1
            print("\nManifold loop")
            print("Iteration Round: {}".format(itt))
            print("Iteration offset: {} (Stops at {})".format(converge, convlimit))

            enthalpies = np.ones(len(self.branches))
            for i in range(len(enthalpies)):
                try:
                    enthalpies[i] = self.branches[i].getStartEnthalpy()
                except:
                    enthalpies[i] = self.refpropm('H','T',self.branches[i].setPointTemp+273.15,'Q',self.branches[i].initialVaporQuality,self.Fluid);

            print("Enthalpy: {} J/kg".format(max(enthalpies)))
            print("Iteration offset: {} (Stops at {})".format(conv_repeat, conv_repeat_limit))
            for ibranch in range(len(self.branches)-1, -1, -1):
               # these will remain unchanged
                self.branches[ibranch].setPointTemp = temperatures[ibranch+1]
                if ibranch >= firstEvapBranch:
                    self.branches[ibranch].finalVaporQualityGuess = vaporQualities[ibranch+1]
                else:
                    self.branches[ibranch].finalUsePressureGuess = True
                    self.branches[ibranch].finalPressureGuess = pressures[ibranch+1]
                # update
                
            print("Running Children")
            p = MyPool(len(self.branches))
            self.branches = p.map(self.worker, self.branches)
            
            # update things that actually change
            for ibranch in range(len(self.branches)-1, -1, -1):
                vaporQualities[ibranch+1] = self.branches[ibranch].getFinalVaporQuality()            
                temperatures[ibranch] = self.branches[ibranch].getInitialTemp()
                pressures[ibranch] = self.branches[ibranch].getStartPressure()

            prev_enthalpies = np.copy(enthalpies)
            for i in range(len(enthalpies)):
                enthalpies[i] = self.branches[i].getStartEnthalpy()
                
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


