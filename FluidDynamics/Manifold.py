class Manifold():
    
    def __init__(self):
        
        self.Name = None
        self.Fluid = None
        self.setPointTemp = None
        self.initialVaporQuality = None
        self.allowedSuperHeatTemp = None 

        self.branches=[]
        self.concatenate = False

