# -*- coding: utf-8 -*-
import sys, os
sys.path.append('PressureDropAndHeatTransfer')
sys.path.append('../REFPROP')
os.environ['RPPREFIX'] = r'C:/Program Files (x86)/REFPROP'  
from Manifold import *

if __name__ == '__main__':

    path = "../XML/HalfBarrel/"
    rootXML = "barrelL1.xml"
    m = Manifold(path, rootXML, 6.8e-3, 0.01, -40., -40.)
    #m.branches[0].setFinalVaporQualityGuess(0.5) # very advanced, don't do this
    m.run()
    m.plot()
