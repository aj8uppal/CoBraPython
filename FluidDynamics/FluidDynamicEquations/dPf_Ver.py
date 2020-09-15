from .Correlations import F_Void_8
from .DarcyWeisbach import *

def dPf_Ver(Fluid, P, H, MFLX, HFLX, Dh, A, Ph, Angle, g, refpropm):
    # RPI = RefPropInterface()
    # refpropm = RPI.refpropm
    dP = -1
    rm = -1
    fluid = Fluid
    q = HFLX

    #refpropm section ADD LATER
    Hl=refpropm('H','P',P,'Q',0,fluid);
    Hv=refpropm('H','P',P,'Q',1,fluid);
    # Density [kg/m3]
    Dl=refpropm('D','P',P,'Q',0,fluid);
    Dv=refpropm('D','P',P,'Q',1,fluid);
    # Viscosity [Pa.s]
    Vl=refpropm('V','P',P,'Q',0,fluid);
    Vv=refpropm('V','P',P,'Q',1,fluid);
    # Thermal Conductivity [W/m.K]
    Kl=refpropm('L','P',P,'Q',0,fluid);
    Kv=refpropm('L','P',P,'Q',1,fluid);
    # Specific Heat [J/kg.K]
    CPl=refpropm('C','P',P,'Q',0,fluid);
    CPv=refpropm('C','P',P,'Q',1,fluid);
    # Surface Tension [N/m]
    ST=refpropm('I','P',P,'Q',0,fluid);

    # Hl = Hv = Dl = Dv = Vl = Vv = Kl = Kv = CPl = CPv = ST = 2 #arbitrary value

    xcrit = 0.6
    ro = 3e-6

    x=refpropm('Q','P',P,'H',H,fluid);
    # x = 2

    #e=F_voidver(x,Dv,g,Dh,Dl,MFLX,ST);
    #above function is not defined, placeholder value
    e = 2

    Tsat=refpropm('T','P',P,'Q',0,Fluid);
    # Tsat = 2
    #[HTC_xt]=F_HTC_xtver(MFLX,Dh,Vl,CPl,Kl);
    HTC_xt = [2]
    #q_ONB=F_q_ONBver(Dv,Tsat,ST,HTC_xt,Hl,Hv,ro);
    q_ONB = 2
    if x == 0:
        #dPf=DarcyWeisbachver(Fluid,P/100,H,MFLX,Dh,A,Ep);
        dPf = DarcyWeisbach(Fluid, P/100, H, MFLX, Dh, A, Ep)
    elif x < xcrit and q > q_ONB:
        # dPf = FriedelCorrelationver(fluid, P, x, MFLX, dh) #arguments not same length (missing A) | kpa/m
        dPf = 2
    else: #x>=xcrit && x<=1
        #mist flow
        dPf = [F_DPmist_45ver(MFLX, Dh, x, Dl, Dv, Vl, Vv)]
        # dPf = FriedelCorrelationver(fluid, P, x, MFLX, dh) #arguments not same length (missing A)
        dPf = 2

    # dPm=F_DPmomentumver(MFLX,q,A,Ph,x,Dl,Dv,Hl,Hv,ST,g,Dh); #not defined
    dPm = 2
    Dtp=Dl*(1-e)+Dv*e;
    # dPsta=F_DPsta(Dtp,g,Angle)      #kpa/m
    dPsta = 2
    dP=(dPf+dPm+dPsta)/100      #bar

    rm=(Dv*e+Dl*(1-e))*A;

    return dP, rm
