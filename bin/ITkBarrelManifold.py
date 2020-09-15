# -*- coding: utf-8 -*-

import sys, os
sys.path.append('../FluidDynamics')
sys.path.append('../REFPROP')
os.environ['RPPREFIX'] = r'C:/Program Files (x86)/REFPROP'  
from Manifold import *

import runinfo

if __name__ == '__main__':

    if len(sys.argv) > 1:
        runinfo.runname = sys.argv[1]
    else:
        runinfo.runname = 'test'    
    
    path = "../XML/HalfBarrel/"
    rootXML = "barrel.xml"
    
    m = Manifold(path, rootXML, 2.*2.304e-3+2.*7.68e-3, 0.0, -1, -40.)
#    m = Manifold(path, rootXML, 2.*2.304e-3, 0.0, -1, -40.)
    m.run()
    m.plot()
