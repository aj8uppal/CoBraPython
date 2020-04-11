from CoolingBranch_v1a import *

class Branch(SingleBranch):
    def __init__(self, Fluid, Tsp, vq0, Tsh, MF, xml, predecessor, successor):
        self.predecessor = predecessor
        self.successor = successor
        super().__init__(Fluid, Tsp, vq0, Tsh, MF, xml)
    def solve(self):
        self.start()
        self.run()
        return {"dP": self.dP[-1]-self.dP[0], "MF": self.massFlow}
    def concatenate(self, other)
        #concatenate branch
        pass
