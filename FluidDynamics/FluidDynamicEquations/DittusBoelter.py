#VERIFIED
# import sys
# sys.path.append('../../REFPROP')
# from refprop import RefPropInterface

def DittusBoelter(Fluid, P, H, MFLX, HFLX, Dh, refpropm):
    '''Input of HFLX is for determining cooling or heating. Ui HFLX is unknown
        then put a -1 for cooling and a +1 for heating.'''
    #Uncomment once refprop is fixed
    # RPI = RefPropInterface()
    # refpropm = RPI.refpropm
    DV=refpropm('V','P',P*1e2,'H',H,Fluid);    # Dynamic Viscosity (Pa*s)
    LAMBDA=refpropm('L','P',P*1e2,'H',H,Fluid)    # W/mK
    CP=refpropm('C','P',P*1e2,'H',H,Fluid) #Specific heat J/kg/K
    DV = 0.2
    LAMBDA = 5
    CP = 2
    PR = DV*CP/LAMBDA #Prandl
    Re = MFLX*Dh/DV

    if HFLX < 0:
        N = 0.33
    else:
        N = 0.4

    NU=0.023*Re**0.8*PR**N;
    HTC=NU*LAMBDA/Dh;
    return HTC
