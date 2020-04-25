# -*- coding: utf-8 -*-
import sys, os
sys.path.append('PressureDropAndHeatTransfer')
sys.path.append('../REFPROP')
os.environ['RPPREFIX'] = r'C:/Program Files (x86)/REFPROP'  
from Manifold import *

if __name__ == '__main__':

    path = "../XML/QuarterShell/"
    rootXML = "endcap.xml"
    m = Manifold(path, rootXML, 17.6e-3/2., 0.01, -35.)
  #  m.branches[0].setFinalVaporQualityGuess(0.5) # very advanced, don't do this
    m.run()
    m.plot()
