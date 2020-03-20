from CoolingBranch_v1a import *
from scipy.interpolate import interp1d

#vary diameter and external heat applied, perhaps length, section of 50cm, one of 1.5m
#1 30W, 1 10W, change dL
#setpointTemp at end must be same
#find the point of the surface that results with the same delta_pressure
#temperature of manifold must be temperature of every branch

#user asks what's the temperature of the inlet or outlet, is temperature everywhere

#at end, pressure is constant, enthalpy is sum

#user passes in setpointtemp, subcoolingtemp, subheattemp

class Manifold:
    def __init__(self, Fluid, Tsp, Tsc, Tsh, MF):
        self.N = 3
        self.Fluid = Fluid
        self.setPointTemp = Tsp #C
        self.subCoolingTemp = Tsc #C
        self.allowedSuperHeatTemp = Tsh #C
        self.massFlow = MF #kg/s
        self.weights = np.full(self.N, self.massFlow/self.N)
        self.dPs = np.zeros_like(self.weights)
        self.epsilon = 1e-3
        self.gamma = 0.00025
        self.totalEnthalpy = -1
        self.iterations = []
    def run(self):
        assert(self.gamma < 0.5)
        self.correction = self.massFlow*1e-1
        self.counts = 0
        # while self.dPs[0] == 0 or abs(max([t - s for s, t in zip(self.dPs, self.dPs[1:])])) > self.epsilon:
            # print(self.epsilon, abs(max([t - s for s, t in zip(self.dPs, self.dPs[1:])])))
        while self.dPs[0] == 0 or self.dPs.std() > self.epsilon:
            self.totalEnthalpy = 0
            for x in range(self.N):
                branch = CoolingBranch_v1a(self.Fluid, self.setPointTemp, self.subCoolingTemp, self.allowedSuperHeatTemp, self.weights[x], "../Manifold.xml", branch=str(x))
                branch.start()
                branch.run(prt=False)
                self.dPs[x] = branch.dP[0]-branch.dP[-1]
                self.totalEnthalpy+=branch.H[-1]
            self.counts+=1
            print("Iteration {}:\n\tmass flows in each branch: {}\n\tdelta pressure in each branch: {}".format(self.counts, self.weights, self.dPs))
            self.refactor()
        else:
            self.present()
    def refactor(self):
        #gradient descent
        # print(self.weights/self.dPs)
        avgDP = self.dPs.mean()
        std = self.dPs.std()
        print("\t\tStandard Deviation of del pressures: {0:.4f}".format(std))
        self.diffs = self.dPs-np.full(self.N, avgDP)
        self.weights-=self.gamma*self.diffs
        # print(self.massFlow, sum(self.weights))
        # breakpoint()
        # self.newweights = list(self.weights)
        # for x in range(self.N-1):
        #     # self.iterations.append([])
        #     ddP = (self.dPs[x+1]-self.dPs[x])/self.dPs[x]
        #     print("\t\tDifference in del pressure between branches {0} and {1}: {2:.4f}%".format(x, x+1, ddP*100))
        #     self.iterations[-1].append(ddP)
        #     # breakpoint()
        #     # ddP/=45
        #     # ddP = sin(max(-pi/2, min(pi/2, ddP)))*0.99
        #     delta = self.weights[x]*ddP*self.gamma
        #     # delta = self.weights[x]*ddP
        #     self.newweights[x]+=delta
        #     self.newweights[x+1]-=delta
        # self.weights = self.newweights
        # ddP = (self.dPs[1]-self.dPs[0])/self.dPs[0]
        # print("\t\tDifference in del pressure: {0:.4f}%".format(ddP*100))
        # self.iterations.append(ddP)
        # delta = self.weights[0]*(ddP*self.alpha)
        # self.weights[0]+=delta
        # self.weights[1]-=delta
    def present(self):
        print("\nConverged.")
        print("End enthalpy: {0:.4f}".format(self.totalEnthalpy))
        for b in range(self.N):
            print("Branch #{}".format(b+1))
            print("\tMass Flow: {}".format(self.weights[b]))
            print("\tDelta Pressure: {}".format(self.dPs[b]))
        print("Initial Mass Flow: {0:.8f}\nEnd Sum Mass Flow: {1:.8f}".format(self.massFlow, sum(self.weights)))
    def plot(self):
        # f = interp1d(range(len(self.iterations)), self.iterations, kind='quadratic')
        # xnew = np.linspace(0, self.iterations[-1], 500)
        pl.plot(range(len(self.iterations)), [i[0] for i in self.iterations])
        pl.xlabel("Iteration #")
        pl.ylabel("Percent difference (dP)")
        pl.show(block=True)

m = Manifold("CO2", -25, 0, 0, 5*1.516e-3)
m.run()
# m.plot()
