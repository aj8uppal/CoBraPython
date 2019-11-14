from Correlations import *
import sys
sys.path.append('../../REFPROP')
from refprop import RefPropInterface

def ThomeCorrelation_EvapHor_CO2(Fluid, P, H, MFLX, HFLX, Dh, A, Ph):
    # function[dP,HTC,x,rm,flowpattern] = ThomeCorrelation_EvapHor_CO2(Fluid,P [bar],H [J/kg],MFLX [kg/m2s], HFLX [ W/m2K], Dh [m], CS-area [m2], periphery [m])
    # Thome correlation.
    fluid=Fluid;
    G=MFLX; #kg/m2s
    d=(4*A/pi)**0.5; # Equivalent diameter.
    p=P;  #bar
    q=HFLX; #w/m2K

    RPI = RefPropInterface()
    refpropm = RPI.refpropm
    #Inlet Vapor Quality
    x=refpropm('Q','P',p*1e2,'H',H,fluid);

    #Enthalpy [J/kg]

    Hl=refpropm('H','P',p*100,'Q',0,fluid);
    Hv=refpropm('H','P',p*100,'Q',1,fluid);

    #Density [kg/m3]
    Dl=refpropm('D','P',p*100,'Q',0,fluid);
    Dv=refpropm('D','P',p*100,'Q',1,fluid);

    #Viscosity [Pa.s]
    Vl=refpropm('V','P',p*100,'Q',0,fluid);
    Vv=refpropm('V','P',p*100,'Q',1,fluid);

    #Thermal Conductivity [W/m.K]
    Kl=refpropm('L','P',p*100,'Q',0,fluid);
    Kv=refpropm('L','P',p*100,'Q',1,fluid);

    #Specific Heat [J/kg.K]
    CPl=refpropm('C','P',p*100,'Q',0,fluid);
    CPv=refpropm('C','P',p*100,'Q',1,fluid);

    #Surface Tension [N/m]
    ST=refpropm('I','P',p*100,'Q',0,fluid);

    # Hl = Hv = Dl = Dv = Vl = Vv = Kl = Kv = CPl = CPv = ST = 2
    # Hv = 2*Hl
    # Dv = 2*Dl
    # x = 4


    #Intermittent to Annular Flow Transition Boundary
    xia=F_xia_18(Dl,Dv,Vl,Vv);
    eia=F_Void_8(G,xia,Dl,Dv,ST);

    #Void fraction
    e=F_Void_8(G,x,Dl,Dv,ST);

    #Get Transition Boundary - Stratified-Wavy to Intermittent and Annular (SW-I/A)
    Gwavy=F_Gwavy_14(G,d,A,x,e,Dl,Dv,ST);
    Gwavy_xia=F_Gwavy_14(G,d,A,xia,eia,Dl,Dv,ST);

    #Get Transition Boundary - Stratified to Stratified-Wavy (S-SW)
    Gstrat=F_Gstrat_17(G,d,A,x,e,Dl,Dv,Vl,Vv,ST);


    #Transition Boundary - Intermittent to Bubbly Flow (I-B)
    Gbub=F_Gbub_26(d,A,x,e,Dl,Dv,Vl);

    #Transition Boundary - Dryout to Mist Flow (D-M)
    Gmist=F_Gmist_24(q,d,x,Dl,Dv,Hl,Hv,ST);

    #Transition Boundary - Annular to Dryout (A-D)
    Gdry=F_Gdry_19(q,d,x,Dl,Dv,Hl,Hv,ST);


    ## flow definition

    #Stratified Flow
    if G<=Gstrat:
        flowpattern='Strat'
    #Slug Flow
    elif G<=Gwavy and x<=xia and G>=Gwavy_xia:
        flowpattern='Slug';

    #Slug+Stratified-Wavy Flow (SLG+SW)
    elif G>=Gstrat and G<=Gwavy and x<=xia and G<Gwavy_xia:
        flowpattern='SlugSW';

    #Stratified-Wavy Flow (SW)
    elif G>=Gstrat and G<=Gwavy and x>=xia and G<=Gdry:
        flowpattern='SW';

    #Intermittent Flow (I)
    elif G>=Gwavy and x<=xia and G<=Gbub:
        flowpattern='Int';

    #Bubbly Flow (B)
    elif G>=Gbub and x<=xia and G<=Gdry:
        flowpattern='Bub';

    #Annular Flow (A)
    elif G>=Gwavy and x>=xia and G<=Gdry:
        flowpattern='Annu';

    #Dry-out (D)
    elif G>=Gstrat and G>=Gdry and G<=Gmist:
        flowpattern='Dry';

    #Mist Flow (M)
    elif G>=Gstrat and G>=Gmist:
        flowpattern='Mist'

    else:
        flowpattern='NotIdentified'
    ##
    if flowpattern == 'Strat':
        DPf=F_DPstrat_52(G,d,x,e,Dl,Dv,Vl,Vv,ST );
        htp=F_HTCstrat( p,G,q,d,A,x,e,CPl,CPv,Dl,Dv,Kl,Kv,Vl,Vv,ST);
    elif flowpattern == 'Slug':
        DPf=F_DPslug_int_35(G,d,x,e,Dl,Dv,Vl,Vv,ST);
        htp=F_HTCannu_int_bub_slug( p,G,q,d,A,x,e,CPl,CPv,Dl,Dv,Kl,Kv,Vl,Vv,ST);
    elif flowpattern == 'SlugSW':
         DPf=F_DPslug_sw_44(G,d,A,x,e,Dl,Dv,Vl,Vv,ST);
         htp=F_HTCslug_sw( p,G,q,d,A,x,e,CPl,CPv,Dl,Dv,Kl,Kv,Vl,Vv,ST);
    elif flowpattern == 'SW':
        DPf=F_DPsw_39(G,d,A,x,e,Dl,Dv,Vl,Vv,ST);
        htp=F_HTCsw( p,G,q,d,A,x,e,CPl,CPv,Dl,Dv,Kl,Kv,Vl,Vv,ST);
    elif flowpattern == 'Int':
        DPf=F_DPslug_int_35(G,d,x,e,Dl,Dv,Vl,Vv,ST);
        htp=F_HTCannu_int_bub_slug( p,G,q,d,A,x,e,CPl,CPv,Dl,Dv,Kl,Kv,Vl,Vv,ST);
    elif flowpattern == 'Bub':
        DPf=F_DPbub_56(G,d,x,e,Dl,Dv,Vl,Vv,ST );
        htp=F_HTCannu_int_bub_slug( p,G,q,d,A,x,e,CPl,CPv,Dl,Dv,Kl,Kv,Vl,Vv,ST);
    elif flowpattern == 'Annu':
        DPf=F_DPannu_29(G,d,x,e,Dl,Dv,Vv,ST);
        htp=F_HTCannu_int_bub_slug( p,G,q,d,A,x,e,CPl,CPv,Dl,Dv,Kl,Kv,Vl,Vv,ST);
    elif flowpattern == 'Dry':
        DPf=F_DPdry_51(G,q,d,A,x,Dl,Dv,Hl,Hv,Vl,Vv,ST);
        htp=F_HTCdry_17( p,G,q,d,A,x,CPl,CPv,Dl,Dv,Hl,Hv,Kl,Kv,Vl,Vv,ST);
    elif flowpattern == 'Mist':
        DPf=F_DPmist_45(G,d,x,Dl,Dv,Vl,Vv);
        htp=F_HTCmist_14( G,d,x,CPv,Dl,Dv,Kv,Vv );
    elif flowpattern == 'NotIdentified':
        DPf=0;
        htp=0;
    end
    DPmom=F_DPmomentum_28(G,q,A,Ph,x,e,Dl,Dv,Hl,Hv,ST);
    #DPmom=0;
    #DPf=0;
    dP=DPf+DPmom;
    ##

    #DP [Pa/m]
    dP=1e-5*dP;

    HTC=htp;
    rm=(Dv*e+Dl*(1-e))*A; #Relative mass kg/m
    return dP, HTC, x, rm, flowpattern
