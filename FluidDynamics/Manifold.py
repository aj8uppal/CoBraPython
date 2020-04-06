from Node import *
from Branch import *

class Manifold():

    def __init__(self):

        self.Name = None
        self.fluid = None
        self.setPointTemp = None
        self.initialVaporQuality = None
        self.allowedSuperHeatTemp = None
        self.initialMassFlow = None
        self.xmls = [None]
        """
                              _________
                     ______v2/   b3    \v3_____
                    /        \_________/       \
                   / b1          b4          b5 \
        v0______v1/                              \v4______v5
            b0    \                              /    b6
                   \                            /
                    \____________b2____________/
        """
        self.branches = {
            0: Branch(self.fluid, self.setPointTemp, vq0, self.allowedSuperHeatTemp, self.massFlow, "xml", predecessor=0, successor=1),
            1: Branch(self.fluid, self.setPointTemp, vq0, self.allowedSuperHeatTemp, self.massFlow, "xml", predecessor=1, successor=2),
            2: Branch(self.fluid, self.setPointTemp, vq0, self.allowedSuperHeatTemp, self.massFlow, "xml", predecessor=1, successor=4),
            3: Branch(self.fluid, self.setPointTemp, vq0, self.allowedSuperHeatTemp, self.massFlow, "xml", predecessor=2, successor=3),
            4: Branch(self.fluid, self.setPointTemp, vq0, self.allowedSuperHeatTemp, self.massFlow, "xml", predecessor=2, successor=3),
            5: Branch(self.fluid, self.setPointTemp, vq0, self.allowedSuperHeatTemp, self.massFlow, "xml", predecessor=3, successor=4),
            6: Branch(self.fluid, self.setPointTemp, vq0, self.allowedSuperHeatTemp, self.massFlow, "xml", predecessor=4, successor=5),
        }
        self.vertices = {
            0: Vertex([self.branches[0]], self.fluid, self.setPointTemp, self.subCoolingTemp, self.allowedSuperHeatTemp, self.massFlow),
            1: Vertex([self.branches[1], self.branches[2]], self.fluid, self.setPointTemp, self.subCoolingTemp, self.allowedSuperHeatTemp, self.massFlow),
            2: Vertex([self.branches[3], self.branches[4]], self.fluid, self.setPointTemp, self.subCoolingTemp, self.allowedSuperHeatTemp, self.massFlow),
            3: Vertex([self.branches[5]], self.fluid, self.setPointTemp, self.subCoolingTemp, self.allowedSuperHeatTemp, self.massFlow),
            4: Vertex([self.branches[6]], self.fluid, self.setPointTemp, self.subCoolingTemp, self.allowedSuperHeatTemp, self.massFlow),
            5: Vertex([], self.fluid, self.setPointTemp, self.subCoolingTemp, self.allowedSuperHeatTemp, self.massFlow),
        }
        # unnecessary?
        # for _, branch in self.branches:
        #     branch.predecessor = self.vertices[branch.predecessor]
        self.root = self.vertices[0]

        # self.branches=[Node(SingleBranch('CO2', -40, 0.01, 0, self.initialMassFlow, self.xmls[0]))]
        self.concatenate = False
    def start(self):
        self.root.initialize()
