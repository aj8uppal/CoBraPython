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
    
    path = "../XML/QuarterShell/"
    rootXML = "endcap.xml"
    
    m = Manifold(path, rootXML, 4*1.219e-3 + 4*3.249e-3 + 8*3.981e-3, 0.0, -1, -40.)
    m.run()
    m.plot()

