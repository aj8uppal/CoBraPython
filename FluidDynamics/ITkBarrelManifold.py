# -*- coding: utf-8 -*-
import sys, os
sys.path.append('PressureDropAndHeatTransfer')
sys.path.append('../REFPROP')
os.environ['RPPREFIX'] = r'C:/Program Files (x86)/REFPROP'  
from Manifold import *

if __name__ == '__main__':

    path = "../XML/HalfBarrel/"
    rootXML = "barrelL0.xml"
    MF = 5e-3
    m = Manifold(path, rootXML, 2e-3, 0.01, -35.)
    m.run()
    m.plot()
