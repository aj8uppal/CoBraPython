# from Node import *
# from SingleBranch import *
from xml.dom import minidom
import numpy as np

class Manifold():

    def __init__(self, base_path, root, initialMassFlow, initialVaporQuality=None, setPointTemp=None, parent=None):

        from SingleBranch import SingleBranch
        self.SingleBranch = SingleBranch
        self.Name = base_path+root
        self.fluid = "CO2"
        self.setPointTemp = setPointTemp
        self.initialVaporQuality = initialVaporQuality or 0.1
        self.allowedSuperHeatTemp = 0
        self.initialMassFlow = initialMassFlow
        # self.xmls = [None]
        self.branches = [] #array of manifolds (branches)
        self.xmls = []
        self.base_path = base_path
        self.root = root
        self.parent = parent
        self.series = None
        self.ST = None
        self.SH = None
        self.dP = 0
        self.initialize()

    def initialize(self):
        # print("Recursively initializing from '{}'...".format(base_path+root))
        self.branches = [(node.getAttribute('type'), node.getAttribute('loc')) for node in minidom.parse(self.base_path+self.root).getElementsByTagName("manifold")]
        self.series = minidom.parse(self.base_path+self.root).getElementsByTagName("components")[0].getAttribute("type") == "series";
        self.setPointTemp = # add XML field which may or may not exist
        self.explore()
    def explore(self):
        for b in range(len(self.branches)):
            branch = self.branches[b]
            if branch[0] == 'manifold':
                self.branches[b] = Manifold(self.base_path, branch[1], self.initialMassFlow if self.series else self.initialMassFlow/len(self.branches), initialVaporQuality=self.initialVaporQuality+(b+1)*(0.5-self.initialVaporQuality)/len(self.branches), setPointTemp=self.setPointTemp, parent=self)
            else:
                self.branches[b] = self.SingleBranch(self.fluid, self.setPointTemp, self.initialVaporQuality if not self.series else self.initialVaporQuality+(b+1)*(0.5-self.initialVaporQuality)/len(self.branches), self.allowedSuperHeatTemp, self.initialMassFlow if self.series else self.initialMassFlow/len(self.branches), self.base_path+branch[1], setPointTemp=self.setPointTemp, parent=self)
    # def solve(self):
    #     for component in self.branches:
    #         if component.__class__ == self.SingleBranch:
    #             component.start(run=True)
    #         else:
    #             component.minimize()
    def updateMassFlow(self, newMF):
        self.initialMassFlow = newMF
        # print(newMF, se
        for branch in self.branches:
            if not self.series:
                branch.updateMassFlow(self.initialMassFlow/len(self.branches))
            else:
                branch.updateMassFlow(self.initialMassFlow)
            # branch.updateMassFlow()
    def minimize(self):
        N = len(self.branches)
        weights = np.full(N, self.initialMassFlow/N)
        print("Total mass flow for {}: {}".format(" ".join(map(str, self.branches)), self.initialMassFlow))
        dPs = np.zeros_like(weights)
        epsilon = 1e-2
        gamma = 0.0001
        totalEnthalpy = -1
        iterations = []
        counts = 0
        # correction = massFlow*1e-1
        # while self.dPs[0] == 0 or abs(max([t - s for s, t in zip(self.dPs, self.dPs[1:])])) > self.epsilon:
            # print(self.epsilon, abs(max([t - s for s, t in zip(self.dPs, self.dPs[1:])])))
        while dPs[0] == 0 or dPs.std() > epsilon:
            totalEnthalpy = 0
            for b in range(len(self.branches)):
                branch = self.branches[b]
                branch.updateMassFlow(weights[b])
                branch.run(run=True)
                dPs[b] = branch.getDP()
                # totalEnthalpy+=branch.getFinalEnthalpy()
            # for x in range(N):
            #     branch = self.SingleBranch(self.Fluid, self.setPointTemp, self.subCoolingTemp, self.allowedSuperHeatTemp, self.weights[x], branches[x].xml)
            #     branch.start()
            #     branch.run(prt=False)
            #     # self.dPs[x] = branch.dP[0]-branch.dP[-1]
            #     branch.getDP()
            #     totalEnthalpy+=branch.getEnthalpy()
            counts+=1
            print("Iteration {}:\n\tmass flows in each branch: {}\n\tdelta pressure in each branch: {}".format(counts, weights, dPs))
            avgDP = dPs.mean()
            std = dPs.std()
            # print("\t\tStandard Deviation of del pressures: {0:.4f}".format(std))
            diffs = dPs-np.full(N, avgDP)
            weights+=gamma*diffs
        return avgDP
    def getInitialTemp(self):
        return self.ST
    def getDP(self):
        return self.dP
    def getFinalVaporQuality(self):
        return self.branches[-1].getFinalVaporQuality()
    def getInitialTemp(self):
        return self.branches[0].getInitialTemp()

    def concat(self):

        vaporQualities = np.linspace(self.initialVaporQuality, 0.5, len(self.branches)+1)
        temperatures = np.ones_like(vaporQualities)*self.setPointTemp

        enthapies = np.ones_like(vaporQualities)*self.refpropm('H','T',self.setPointTemp+273.15,'Q',self.initialVaporQuality,self.Fluid)

        itt=0;
        converge=10000;
        convlimit=40;
        conv_repeat=0;
        conv_repeat_limit=3;
        ittstop = 400

        while (abs(converge)>convlimit or conv_repeat<conv_repeat_limit+1) and itt<ittstop:
            itt+=1

            for ibranch in range(len(branches)-1, -1, -1):
                # these will remain unchanged
                self.branches[ibranch].setPointTemp = temperatures[ibranch+1]
                self.branches[ibranch].initialVaporQuality = vaporQualities[ibranch]
                # this will change
                self.branches[ibranch].finalVaporQualityGuess = vaporQualities[ibranch+1]
                # update
                self.branches[ibranch].run(run = True)
                # update things that actually change
                vaporQualities[ibranch+1] = self.branches[ibranch].getFinalVaporQuality()            
                temperatures[ibranch] = self.branches[ibranch].getInitialTemp()

            prev_enthalpies = np.copy(entalpies)
            for i,enthalpy in enumerate(enthalpies):
                enthalpy = self.refpropm('H','T',temperatures[i]+273.15,'Q',vaporQualities[i],self.Fluid)

            diff_enthalpies = np.abs(enthalpies-prev_enthalpies)
            converge = np.amax(diff_enthalpies) #use enthalpy to converge
            if abs(converge) < convlimit:
                conv_repeat+=1
            else:
                conv_repeat = 0

        return vaporQualities[-1]
    
    # def __repr__(self):
        # return self.Name
    def __str__(self):
        return self.root
    # def prettyprint(self):
    #     if self.branches == None:
    #         return
    #     print(self.Name)
    #     if self.series:
    #         print("\t".join(map(str, self.branches)))
    #     else:
    #         print("\n".join(map(str, self.branches)))
        # print("\t".join([b if b.__class__ == self.SingleBranch else b.prettyprint() for b in self.branches]))
    def prettyprint(self, i=0, var=None):
        print(i*"\t|"+self.Name, end="    ")
        print(self.initialMassFlow if var == "IMF" else "")
        for b in self.branches:
            if b.__class__ == Manifold:
                b.prettyprint(i=i+1, var=var)
            else:
                print((i+1)*"\t|"+str(b), end="    ")
                print(b.massFlow if var == "IMF" else "")
    def run(self, ST=None, run=None):
        ST = ST or self.setPointTemp
        if self.series == True:
            print("Concatenating {}".format(" ".join(map(str, self.branches))))
            return self.concat()
            # print(" + ".join(map(str, self.branches)))
        else:
            print("Minimizing {}".format(" ".join(map(str, self.branches))))
            # print("DP Minimization on {}".format(" ".join(map(str, self.branches))))
            # if all([type(branch) == self.SingleBranch for branch in self.branches]):
            #     return self.minimize(self.branches, 1)
            #     #self.solve(self.branches)
            # else:
            return self.minimize()
            # return self.minimize(*[branch.run() if type(branch) == Manifold else branch for branch in self.branches])
            # print("Minimization on "+", ".join(map(str, self.branches)))
            # return 3
            # print(" / ".join(map(str, self.branches)))
        # for branch in self.branches:
        #     branch.run()
    # def run(self):
    #     for branch in branches:
    #
    #         endpoint = branch
    #         if all(branch.successor)
    # def start(self):
    #     self.root.initialize()



# b3b4 = Manifold()

path = "../XML/Manifoldv1/"
rootXML = "manifold0.xml"
MF = 5e-3
m = Manifold(path, rootXML, MF)
# b = m.initialize(path, "manifold0.xml")
# root = m.restructure(m.root, m.end)
# print(root.run())
# m.prettyprint()
# n = Manifold(path, root="branch1.xml")
# n.prettyprint()
# o = Manifold(path, root="branch3.xml")
# o.prettyprint()

"""
                    Fig. 1
                  __________
           ______/    b3    \_____
          /      \__________/     \
         / b1         b4        b5 \
  ______/                           \______
    b0  \                           /  b6
         \                         /
          \___________b2__________/
            Run b5, guess initial vapor quality, gives vapor quality that we guessed and Temperature


                    Fig. 2

____b0__________b2__________b4_________b6_____
-->       |           |           |           |
          |           |           |           |
          |           |           |           |
        b1|         b3|         b5|         b6|
          |           |           |           |
          |           |           |           |
          |           |           |           |
____b9____|____b8_____|_____b7____|____b6_____|
<--
"""
