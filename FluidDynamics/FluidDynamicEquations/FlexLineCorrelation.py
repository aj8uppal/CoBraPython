#from Correlations import F_Void_8
import sys
# sys.path.append('../../REFPROP')
# from refprop import RefPropInterface

#VERIFIED

def FlexLineCorrelation(Fluid, MFLX):
    # RPI = RefPropInterface()
    # refpropm = RPI.refpropm
    #Fluid properties: (arbitrary values until refprop is fixed)

    # This is what is measured by the strips
    dP = (99.47-MFLX)/30.45
    rm = 0
    return dP, rm
