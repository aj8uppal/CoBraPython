# from Node import *
# from SingleBranch import *
from xml.dom import minidom

class Manifold():

    def __init__(self, branches=[], concatenate=False):

        from SingleBranch import SingleBranch
        self.SingleBranch = SingleBranch
        self.Name = None
        self.fluid = "CO2"
        self.setPointTemp = -40
        self.initialVaporQuality = 0.01
        self.allowedSuperHeatTemp = 0
        self.initialMassFlow = None
        self.xmls = [None]
        self.branches = branches #array of manifolds (branches)
        self.xmls = []
        self.concatenate = concatenate

    def initialize(self, base_path, root):
        self.base_path = base_path
        self.root = root
        print("Recursively initializing from '{}'...".format(base_path+root))

        self.explore(root)
        # self.root = [b for b in self.branches if b.xml == base_path+root][0]
        self.root = [b for b in self.branches if not len(b.predecessors)][0]
        self.end = [b for b in self.branches if not b.successors[0]][0]
        #self.initialize_loops(self.root)
    def explore(self, rootXML, parent=None):
        if rootXML == "None":
            return None
        rootXML = self.base_path+rootXML
        if rootXML in self.xmls:
            branch = [b for b in self.branches if b.xml == rootXML][0]
            branch.predecessors.append(parent)
            # print("\t".join([str(branch), str(list(map(str, branch.predecessors))), str(list(map(str, branch.successors)))]))
            return branch
        self.xmls.append(rootXML)
        successors = minidom.parse(rootXML).getElementsByTagName("data")[0].getAttribute("successors").split(",")
        branch = self.SingleBranch(self.fluid, self.setPointTemp, self.initialVaporQuality, self.allowedSuperHeatTemp, self.initialMassFlow, rootXML, predecessors=[parent])
        if not branch.predecessors[0]:
            branch.predecessors = []
        branch.successors = [self.explore(child_xml, branch) for child_xml in successors]
        # print("\t".join([str(branch), str(list(map(str, branch.predecessors))), str(list(map(str, branch.successors)))]))
        self.branches.append(branch)
        return branch
    def restructure(self, branch, end):
        if branch is None:
            return
        preds = len(branch.predecessors)
        succs = len(branch.successors)
        #CASE 1: Series-Parallel
        if not preds: #root
            return Manifold(branches=[branch, Manifold(branches=[self.restructure(successor, end) for successor in branch.successors], concatenate=False), end], concatenate=True)
        #CASE 2: Parallel
        if succs > 1:
            if len(branch.successors[0].successors) == 1 and succs == len(branch.successors[0].successors[0].predecessors): #simplifiable
                return Manifold(branches=[branch, Manifold(branches=branch.successors, concatenate=False), branch.successors[0].successors[0]], concatenate=True)
        #CASE 3: Series
        else:
            return branch
    def initialize_loops(self, b):
        if b is None:
            return
        if b.successors[0] and b.successors[0].successors[0] and len(b.successors) == len(b.successors[0].successors[0].predecessors):
            print(", ".join(map(str, b.successors))+" forms a simplifiable loop")
        # if len(b.successors) == 1 and b.successors[0].predecessors == [branch for branch in self.branches if b.successors[0] in branch.successors]:
        #     print("Simplifiable loop, ")
        for i in b.successors:
            self.initialize_loops(i)
        # if len(b.predecessors) == 1:
        #     if len(b.predecessor[0].successors) == 1):
        #         print("{} is a single branch".format(b.xml))
        #     elif b.
        # if (len(b.predecessors) == 1 and len(b.predecessor[0].successors) > 1)
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
m = Manifold()
m.initialize(path, root="branch0.xml")
root = m.restructure(m.root, m.end)
print(root.run())
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
