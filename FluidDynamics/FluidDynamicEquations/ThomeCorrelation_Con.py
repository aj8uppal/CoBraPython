from .Correlations import *
from numpy import isnan

def ThomeCorrelation_Con(Fluid, P, H, MFLX, HFLX, Dh, A, Ph, refpropm):
    # RPI = RefPropInterface()
    # refpropm = RPI.refpropm
    fluid = Fluid
    G=MFLX #kg/m2s
    d=Dh #Equivalent diameter.
    p=P*1e2
    q=abs(HFLX) #w/m2K
    #Commented until refprop is fixed, arbitrary value for now
    # Inlet Vapor Quality
    x=abs(refpropm('Q','P',p*100,'H',H,fluid));

    #Enthalpy [J/kg]

    Hl=refpropm('H','P',p,'Q',0,fluid);
    Hv=refpropm('H','P',p,'Q',1,fluid);

    #Density [kg/m3]
    Dl=refpropm('D','P',p,'Q',0,fluid);
    Dv=refpropm('D','P',p,'Q',1,fluid);

    #Viscosity [Pa.s]
    Vl=refpropm('V','P',p,'Q',0,fluid);
    Vv=refpropm('V','P',p,'Q',1,fluid);

    #Thermal Conductivity [W/m.K]
    Kl=refpropm('L','P',p,'Q',0,fluid);
    Kv=refpropm('L','P',p,'Q',1,fluid);

    #Specific Heat [J/kg.K]
    CPl=refpropm('C','P',p,'Q',0,fluid);
    CPv=refpropm('C','P',p,'Q',1,fluid);

    #Surface Tension [N/m]
    ST=refpropm('I','P',p,'Q',0,fluid);

    # x = Hl = Hv = Dl = Dv = Vl = Vv = Kl = Kv = CPl = CPv = ST = 2
    # Hv = 2*Hl
    #Intermittent to Annular Flow Transition Boundary
    xia= F_xia_18con(Dl,Dv,Vl,Vv);
    eia=F_Void_8(G,xia,Dl,Dv,ST);

    #Void fraction
    e=F_Void_8(G,x,Dl,Dv,ST);

    #Get Transition Boundary - Stratified to Stratified-Wavy (S-SW)
    Gstrat=F_Gstrat_17con(d,A,x,e,Dl,Dv,Vl);


    #Get Transition Boundary - Stratified-Wavy to Intermittent and Annular (SW-I/A)
    # print(x, Dl, Dv)
    Gwavy=F_Gwavy_16con(G,d,A,x,e,Dl,Dv,ST);

    #Transition Boundary - Intermittent to Bubbly Flow (I-B)
    Gbub=F_Gbub_22con(d,A,x,e,Dl,Dv,Vl);

    #Transition Boundary - Annular and intermittent flow to Mist Flow (AI-M)
    Gmist=F_Gmist_19con(G,d,A,x,e,Dl,Dv,ST);
    # print(Gmist, x, xia, G, Gbub)
    #flow definition

    # print(G, Gmist, x, xia)
    if G <= Gstrat:
        #fully stratified flow
        flowpattern = 'Strat'
        FilmAngle = F_StratAngle_13(e)
    elif G >= Gstrat and G <= Gwavy:
        flowpattern = 'SW'
        FilmAngle = F_FilmAngle_8130con(G, d, A, x, e, Dl, Dv, Vl, ST)
    elif G >= Gwavy and G <= Gmist.real and x <= xia and G <= Gbub.real:
        flowpattern='Int';
        FilmAngle = 0
    elif G <= Gmist.real and x <= xia and G > Gbub.real:
        flowpattern = 'Bub'
        FilmAngle = 0
    elif G>=Gwavy and G<=Gmist.real and x>xia:
        flowpattern = 'Annu'
        FilmAngle = 0
    elif G > Gmist.real and G >= Gstrat:
        flowpattern = 'Mist'
        FilmAngle = 0
    else:
        flowpattern = 'NotIdentified'
    # print(flowpattern)
    if flowpattern == 'Strat':
        DPf=F_DPstrat_52(G,d,x,e,Dl,Dv,Vl,Vv,ST );
        htp=F_HTC_8123con( G,q,d,A,x,e,FilmAngle,Dl,Dv,CPl,Kl,Hl,Hv,Vl,ST)
    elif flowpattern == 'SW':
        DPf=F_DPsw_39(G,d,A,x,e,Dl,Dv,Vl,Vv,ST);
        htp=F_HTC_8123con( G,q,d,A,x,e,FilmAngle,Dl,Dv,CPl,Kl,Hl,Hv,Vl,ST);
    elif flowpattern == 'Int':
        DPf=F_DPslug_int_35(G,d,x,e,Dl,Dv,Vl,Vv,ST);
        htp=F_HTC_8123con( G,q,d,A,x,e,FilmAngle,Dl,Dv,CPl,Kl,Hl,Hv,Vl,ST);
        # print(G, q, d, A, x, e, FilmAngle, Dl, Dv, CPl, Kl, Hl, Hv, Vl, ST, p, H);
    elif flowpattern == 'Annu':
        DPf=F_DPannu_29(G,d,x,e,Dl,Dv,Vv,ST);
        htp=F_HTC_8123con( G,q,d,A,x,e,FilmAngle,Dl,Dv,CPl,Kl,Hl,Hv,Vl,ST);
    elif flowpattern == 'Mist':
        DPf=F_DPmist_45(G,d,x,Dl,Dv,Vl,Vv);
        htp=F_HTC_8123con( G,q,d,A,x,e,FilmAngle,Dl,Dv,CPl,Kl,Hl,Hv,Vl,ST);
        # if(isnan(htp)):
        #     breakpoint()
    elif flowpattern == 'Bub':
        DPf=F_DPbub_56(G,d,x,e,Dl,Dv,Vl,Vv,ST );
        htp=F_HTC_8123con( G,q,d,A,x,e,FilmAngle,Dl,Dv,CPl,Kl,Hl,Hv,Vl,ST);
    elif flowpattern == 'NotIdentified':
        DPf=0
        htp=0

    DPmom=F_DPmomentum_28(G,q,A,Ph,x,e,Dl,Dv,Hl,Hv,ST);
    dP=DPf+DPmom;

    dP=1e-5*dP #bar
    HTC=htp.real;

    rm=(Dv*e+Dl*(1-e))*A #Relative mass kg/m
    return dP, HTC, x, rm, flowpattern, xia, Gwavy, None, Gstrat, Gbub, Gmist, None
