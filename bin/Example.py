# -*- coding: utf-8 -*-
import sys, os
sys.path.append('PressureDropAndHeatTransfer')
sys.path.append('../REFPROP')
os.environ['RPPREFIX'] = r'C:/Program Files (x86)/REFPROP'  
from Manifold import *

if __name__ == '__main__':

    path = "../XML/Manifoldv1/"
    rootXML = "manifold0.xml"
    MF = 5e-3
    m = Manifold(path, rootXML, 2.02e-3, 0.01, -35.)
    m.run()
    m.plot()
"""
                    Fig. 1
                  __________
           ______/    b3    \_____
          /      \__________/     \
         / b1         b4        b5 \
  ______/                           \______
    b0  \                           /  b6
         \                         /
          \___________b2__________/
            Run b5, guess initial vapor quality, gives vapor quality that we guessed and Temperature


                    Fig. 2

____b0__________b2__________b4_________b6_____
-->       |           |           |           |
          |           |           |           |
          |           |           |           |
        b1|         b3|         b5|         b6|
          |           |           |           |
          |           |           |           |
          |           |           |           |
____b9____|____b8_____|_____b7____|____b6_____|
<--
"""
