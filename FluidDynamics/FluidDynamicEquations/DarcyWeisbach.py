# import sys
# sys.path.append('../../REFPROP')
# from refprop import RefPropInterface
from ColebrookEquation import *
from BlasiusCorrelation import *
from time import time

#VERIFIED

def DarcyWeisbach(Fluid, P, H, MFLX, Dh, A, Ep, refpropm):
    '''Darcy Weissbach formula for single phase pressure drop
    Inputs: Fluid, P (bar), H (kJ/kg), MFLX (kg/m2), Dh (m)
    Output: dP (bar)'''
    #Uncomment these when refprop is fixed
    # RPI = RefPropInterface()
    # refpropm = RPI.refpropm
    # start = time()
    # print(P.real*1e2, H, Fluid)
    
    #section takes a lot more time if P is large...
    RHO=refpropm('D','P',P*1e2,'H',H,Fluid) #kg/m^3
    DV=refpropm('V','P',P*1e2,'H',H,Fluid) #dynamic viscosity
    # print(time()-start)
    RHO = 5
    DV = 0.2

    Re = MFLX*Dh/DV
    if Re<4000:
        (Fd,State)=BlasiusCorrelation(Re)
    else:
        Fd=ColebrookEquation(Dh,Ep,Re)
        State='Tur'

    dP=Fd*MFLX**2/(2*Dh*RHO)/1e5 #bar/m
    rm=RHO*A #Relative mass kg/m
    return dP, rm, State
