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
        self.N = 2
        self.Fluid = Fluid
        self.setPointTemp = Tsp #C
        self.subCoolingTemp = Tsc #C
        self.allowedSuperHeatTemp = Tsh #C
        self.massFlow = MF #kg/s
        self.weights = np.full(self.N, self.massFlow/self.N)
        self.dPs = np.zeros_like(self.weights)
        self.epsilon = 1e-6
        self.alpha = 0.25
        self.div = 1000
        self.totalEnthalpy = -1
        self.iterations = []
    def run(self):
        assert(self.alpha < 0.5)
        self.correction = self.massFlow*1e-1
        self.counts = 0
        while self.dPs[0] == 0 or abs(max([t - s for s, t in zip(self.dPs, self.dPs[1:])])) > self.epsilon:
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
        ddP = (self.dPs[1]-self.dPs[0])/self.dPs[0]
        print("\t\tDifference in del pressure: {0:.4f}%".format(ddP*100))
        self.iterations.append(ddP)
        self.weights[0]*=(1+(ddP*self.alpha))
        self.weights[1]/=(1+(ddP*self.alpha))
    def present(self):
        print("\nConverged.")
        print("End enthalpy: {0:.4f}".format(self.totalEnthalpy))
        for b in range(self.N):
            print("Branch #{}".format(b+1))
            print("\tMass Flow: {}".format(self.weights[b]))
            print("\tDelta Pressure: {}".format(self.dPs[b]))
    def plot(self):
        f = interp1d(range(len(self.iterations)), self.iterations, kind='cubic')
        pl.plot(range(len(self.iterations)), f(range(len(self.iterations))))
        pl.xlabel("Iteration #")
        pl.ylabel("Percent difference (dP)")
        pl.show(block=True)

m = Manifold("CO2", -25, 0, 0, 5*1.516e-3)
m.run()
m.plot()
