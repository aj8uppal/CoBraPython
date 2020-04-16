# from Node import *
# from SingleBranch import *
from xml.dom import minidom

class Manifold():

    def __init__(self, base_path, root, parent=None):

        from SingleBranch import SingleBranch
        self.SingleBranch = SingleBranch
        self.Name = None
        self.fluid = "CO2"
        self.setPointTemp = -40
        self.initialVaporQuality = 0.01
        self.allowedSuperHeatTemp = 0
        self.initialMassFlow = None
        # self.xmls = [None]
        self.branches = [] #array of manifolds (branches)
        self.xmls = []
        self.base_path = base_path
        self.root = root
        self.parent = parent
        self.initialize()

    def initialize(self):
        # print("Recursively initializing from '{}'...".format(base_path+root))
        self.branches = [(node.getAttribute('type'), node.getAttribute('loc')) for node in minidom.parse(self.base_path+self.root).getElementsByTagName("manifold")]
        self.series = minidom.parse(self.base_path+self.root).getElementsByTagName("components")[0].getAttribute("type") == "series";
        self.explore()
    def explore(self):
        for b in range(len(self.branches)):
            branch = self.branches[b]
            if branch[0] == 'manifold':
                self.branches[b] = Manifold(self.base_path, branch[1], parent=self)
            else:
                self.branches[b] = self.SingleBranch(self.fluid, self.setPointTemp, self.initialVaporQuality, self.allowedSuperHeatTemp, self.initialMassFlow, self.base_path+branch[1], parent=self)
    def solve(self):
        for component in self.branches:
            if component.__class__ == self.SingleBranch:
                component.start(run=True)
            else:
                component.minimize()
    def minimize(self, branches, massFlow):
        return 1
        N = len(branches)
        weights = np.full(N, massFlow/N)
        dPs = np.zeros_like(weights)
        epsilon = 1e-3
        gamma = 0.0001
        totalEnthalpy = -1
        iterations = []
        correction = massFlow*1e-1
        counts = 0
        # while self.dPs[0] == 0 or abs(max([t - s for s, t in zip(self.dPs, self.dPs[1:])])) > self.epsilon:
            # print(self.epsilon, abs(max([t - s for s, t in zip(self.dPs, self.dPs[1:])])))
        while dPs[0] == 0 or dPs.std() > epsilon:
            totalEnthalpy = 0
            for x in range(N):
                branch = self.SingleBranch(self.Fluid, self.setPointTemp, self.subCoolingTemp, self.allowedSuperHeatTemp, self.weights[x], branches[x].xml)
                branch.start()
                branch.run(prt=False)
                # self.dPs[x] = branch.dP[0]-branch.dP[-1]
                branch.getDP()
                totalEnthalpy+=branch.getEnthalpy()
            counts+=1
            # print("Iteration {}:\n\tmass flows in each branch: {}\n\tdelta pressure in each branch: {}".format(self.counts, self.weights, self.dPs))
            avgDP = Ps.mean()
            std = dPs.std()
            # print("\t\tStandard Deviation of del pressures: {0:.4f}".format(std))
            diffs = dPs-np.full(N, avgDP)
            weights-=gamma*diffs
        return weights, totalEnthalpy
    def concat(self, *args):
        return 2
        # return sum([*args])
    def solve(self, branches):
        pass
    def prettyprint(self):
        for branch in self.branches:
            # print("{} splits from\n".format(branch.xml)+"\n".join(["\t"+p.xml for p in branch.predecessors if p]))
            print("{} splits into\n".format(branch.xml)+"\n".join(["\t"+s.xml for s in branch.successors if s]))
    def run(self):
        if self.concatenate == True:
            print("Concatenation on {}".format(" ".join(map(str, self.branches))))
            return self.concat(*[branch.run() for branch in self.branches])
            # print(" + ".join(map(str, self.branches)))
        else:
            print("DP Minimization on {}".format(" ".join(map(str, self.branches))))
            if all([type(branch) == self.SingleBranch for branch in self.branches]):
                return self.minimize(self.branches, 1)
                #self.solve(self.branches)
            else:
                return self.minimize(*[branch.run() if type(branch) == Manifold else branch for branch in self.branches])
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
m = Manifold(path, "manifold0.xml")
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
