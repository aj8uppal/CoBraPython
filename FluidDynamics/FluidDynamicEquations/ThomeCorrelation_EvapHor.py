from Correlations import *
# import sys
# sys.path.append('../../REFPROP')
# from refprop import RefPropInterface

def ThomeCorrelation_EvapHor(Fluid, P, H, MFLX, HFLX, Dh, A, Ph, refpropm):
    # Wojtan-Ursenbacher-Thome model_updated based on model of
    # Kattan-Thome-Favrat Map_Databook  Chapter 12 flow pattern map in
    # horizontal plain tube
    # RPI = RefPropInterface()
    # refpropm = RPI.refpropm
    fluid=Fluid;
    G=MFLX; #kg/m2s
    # d=(4*A/pi)^0.5; # Equivalent diameter.
    d=Dh;
    p=P;  #bar
    q=HFLX; #w/m2K
    g=9.81;

    #Inlet Vapor Quality
    x=refpropm('Q','P',round(p*1e2),'H',round(H),fluid);

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


    Hl = Hv = Dl = Dv = Vl = Vv = Kl = Kv = CPl = CPv = ST = 2
    Hv = 2*Hl

    #Intermittent to Annular Flow Transition Boundary
    xia= F_xia_18con(Dl,Dv,Vl,Vv);  # it is the same with that of condensation
    eia=F_Void_8(G,xia,Dl,Dv,ST);
    Gwavy_xia=F_Gwavy_14(G,d,A,xia,eia,Dl,Dv,ST);

    #Void fraction
    e=F_Void_8(G,x,Dl,Dv,ST);

    #Get Transition Boundary - Stratified to Stratified-Wavy (S-SW)
    Gstrat= F_Gstrat_12417evap(G,d,A,x,e,Dl,Dv,Vl,Vv,ST);
    # fprintf('Gstrat is #f\n:',Gstrat);

    #Get Transition Boundary - Stratified-Wavy to Intermittent and Annular (SW-I/A)
    Gwavy=F_Gwavy_14(G,d,A,x,e,Dl,Dv,ST);

    #Transition Boundary - Intermittent to Bubbly Flow (I-B)
    # Gbub=F_Gbub_26(d,A,x,e,Dl,Dv,Vl);   # The formular is the same with that of CO2 evaporation
    # bubbly not mentioned in databook

    #Transition Boundary - Annular to Dryout (A-D)
    Gdry=F_Gdry_12429evap(q,d,x,Dl,Dv,Hl,Hv,ST);

    #Transition Boundary - Dryout to Mist Flow (D-M)
    Gmist=F_Gmist_12430evap(q,d,x,Dl,Dv,Hl,Hv,ST);


    ## flow definition

    #Stratified Flow
    if G<=Gstrat:
        flowpattern='Strat'
    #Slug Flow
    elif G<=Gwavy and x<=xia and G>=Gwavy_xia:
        flowpattern='Slug'
    #Slug+Stratified-Wavy Flow (SLG+SW)
    elif G>Gstrat and G<Gwavy and x<xia and G<Gwavy_xia:
        flowpattern='SlugSW'
    #Stratified-Wavy Flow (SW)
    elif G>Gstrat and G<=Gwavy and x>=xia and G<=Gdry:
        flowpattern='SW'
    #Intermittent Flow (I)
    elif G>=Gwavy and x<=xia: #and G<=Gbub
        flowpattern='Int'
    #Bubbly Flow (B)
    #elseif G>=Gbub and x<=xia and G<=Gdry
    #flowpattern='Bub';
    #Annular Flow (A)
    elif G>=Gwavy and x>xia and G<=Gmist and G<=Gdry:
        flowpattern='Annu'

    #Dry-out (D)
    elif G>=Gstrat and G>=Gdry and G<=Gmist:
        flowpattern='Dry'
    #Mist Flow (M)
    elif G>=Gstrat and G>Gmist:
        flowpattern='Mist'
    else:
        flowpattern='NotIdentified';
    ##
    if flowpattern == 'Strat':
        #DPf=F_DPsta_evap(G,d,x,e,Dl,Dv,Vv,Vl,ST,q,Hl,Hv);  #done
        htp=F_HTCstrat_evap( fluid,G,q,d,A,x,e,CPl,CPv,Kl,Kv,Vl,Vv); #done
    elif flowpattern == 'SW':
        #DPf=F_DPsw_evap(G,d,x,e,Dl,Dv,Vv,Vl,A,ST);#done
        htp=F_HTC_sw_evap(fluid,G,q,d,A,x,e,Dl,Dv,CPl,CPv,Hl,Hv,Kl,Kv,Vl,Vv,ST); #done
    elif flowpattern == 'SlugSW':
         #DPf= F_DPslug_sw_13256_evap(G,d,x,e,eia,Dl,Dv,Vv,Vl,A,ST); #done
         htp=F_HTCslug_sw_evap(fluid,G,q,d,A,x,e,Dl,Dv,CPl,CPv,Kl,Kv,Hl,Hv,Vl,Vv,ST); #done
    elif flowpattern == 'Slug':
        #DPf=F_DPslug_int_13249_evap(G,d,x,e,Dl,Dv,Vl,Vv,ST,A);#done
        htp=F_HTCannu_int_bub_slug_evap(fluid,G,q,d,A,x,e,CPl,CPv,Kl,Kv,Vl,Vv); #done
    elif flowpattern == 'Int':
       # DPf=F_DPslug_int_13249_evap(G,d,x,e,Dl,Dv,Vl,Vv,ST,A);#done
        htp=F_HTCannu_int_bub_slug_evap(fluid,G,q,d,A,x,e,CPl,CPv,Kl,Kv,Vl,Vv); #done
    # elseif strcmp(flowpattern,'Bub')==1  #It is not addressed
    #     DPf=F_DPbub_56(G,d,x,e,Dl,Dv,Vl,Vv,ST );
    #        htp=F_HTCannu_int_bub_slug_evap(fluid,G,q,d,A,x,e,CPl,Cpv,Kl,Kv,Vl,Vv); #done
    elif flowpattern == 'Annu':
       #DPf=F_DPannu_13244_evap(G,d,x,e,Dl,Dv,Vv,Vl,ST,A); #done
        htp=F_HTCannu_int_bub_slug_evap(fluid,G,q,d,A,x,e,CPl,CPv,Kl,Kv,Vl,Vv);
    elif flowpattern == 'Dry':
        #DPf=F_DPdry_13258_evap(G,d,x,e,Dl,Dv,Vv,Vl,Hl,Hv,ST,q,g,A); #done
        htp=F_HTCdry_18711evap( fluid,G,q,d,A,x,CPl,CPv,Dl,Dv,Kl,Kv,Vl,Vv,Hl,Hv,ST); #done
    elif flowpattern == 'Mist':
        #DPf=F_DPmist_evap(G,d,x,Dl,Dv,Vl,Vv);#done
        htp=F_HTCmist_18710evap( G,d,x,CPv,Dl,Dv,Kv,Vv ); #done
    elif flowpattern == 'NotIdentified':
        DPf=0;
        htp=0;

    DPf=0;

    DPmom=F_DPmomentum_28(G,q,A,Ph,x,e,Dl,Dv,Hl,Hv,ST);
    dP=DPf+DPmom;
    ##
    dP=dP*(1e-5); #bar
    HTC=htp;
    rm=(Dv*e+Dl*(1-e))*A; #Relative mass kg/m
    return dP,HTC,x,rm,flowpattern
