from cmath import log, sqrt, log10
from numpy import pi, cos, exp, sin
import sys
sys.path.append('..')
from BlasiusCorrelation import *

def F_Void_8(G, x, Dl, Dv, ST):
    g = 9.81
    # Void Fraction (Rouhani-Axelsson drift flux model)
    #e=(x/Dv)*((1+0.12*(1-x))*((x/Dv)+((1-x)/Dl))+((1.18*(1-x)*(g*ST*(Dl-Dv))**(1/4))/(G*Dl**(1/2))))**(-1)

    # Homogeneous Void Fraction
    eH=1/(1+((1-x)/x)*(Dv/Dl));

    #Rouhani Void Fraction
    er=(x/Dv)*((1+0.12*(1-x))*((x/Dv)+((1-x)/Dl))+((1.18*(1-x)*(g*ST*(Dl-Dv))**(1/4))/(G*Dl**(1/2))))**(-1)

    #Logarithmic Mean Void Fraction
    # e = (eH-er)/log(eH/er)
    # e = er
    return er

def F_xia_18(Dl, Dv, Vl, Vv):
    #Intermittent to Annular Flow Transition Boundary
    xia=(1.8**(1/0.875)*(Dv/Dl)**(-1/1.75)*(Vl/Vv)**(-1/7)+1)**(-1);
    return xia

def F_xia_18con(Dl, Dv, Vl, Vv):
    xia=(0.2914*(Dv/Dl)**(-1/1.75)*(Vl/Vv)**(-1/7)+1)**(-1) #0.2914 replaced 1.8**(1/0.875)
    return xia

def F_Ald_9(d, A, e):
    Ald = A*(1-e)/d**2
    return Ald

def F_Avd_10(d, A, e):
    Avd = A*e/d**2
    return Avd

def F_Gstrat_17con(d, A, x, e, Dl, Dv, Vl):
    g = 9.81

    #Dimensionless Cross Sectional Area Occupied by Vapor Phase
    Avd=F_Avd_10( d,A,e );
    #Dimensionless Cross Sectional Area Occupied by Liquid Phase
    Ald=F_Ald_9( d,A,e );

    #Transition Boundary - Stratified to Stratified-Wavy (S-SW)
    # if x<xia;
    # Gstrat=((226.3**2*Ald_xia*Avd_xia**2*Dv*(Dl-Dv)*Vl*g)/(xia**2*(1-xia)*pi**3))**(1/3);
    # else
    Gstrat=((226.3**2*Ald*Avd**2*Dv*(Dl-Dv)*Vl*g)/(x**2*(1-x)*pi**3))**(1/3)+20*x;   #20*x is new
    return Gstrat

def F_StratAngle_13(e):
    Ang_strat=2*pi-2*(pi*(1-e)+((3*pi)/2)**(1/3)*(1-2*(1-e)+(1-e)**(1/3)-e**(1/3))-(1/200)*(1-e)*e*(1-2*(1-e))*(1+4*((1-e)**2+e**2)))
    return Ang_strat

def F_hld_11(e):
    StratAngle=F_StratAngle_13(e)
    hld=0.5*(1-cos((2*pi-StratAngle)/2))
    return hld

def F_Gwavy_16con(G, d, A, x, e, Dl, Dv, ST):
    #Kattan-Thome-Favrat criterion
    g = 9.81
    #G=1; # G as fake input to weber and Froude, will be divided away.

    #Dimensionless Cross Sectional Area Occupied by Vapor Phase
    Avd=F_Avd_10( d,A,e );

    #Dimensionless Vertical Height of the Liquid
    hld = F_hld_11(e);

    #Liquid Froude Number
    #Frl=F_FroudeLiquid_15( G,d,Dl );

    #Liquid Weber Number for Flow Pattern
    #Wel=F_Weber_16_21( G,d,Dl,ST);

    Wel_Frl=g*d**2*Dl/ST #(We/Fr)l
    Gwavy=(((16*Avd**3*g*d*Dl*Dv)/(x**2*pi**2*(1-(2*hld-1)**2)**(1/2)))*((pi**2/(25*hld**2))*(Wel_Frl)**(-1.023)+1))**(1/2)+50-75*exp(-(x**2-0.97)**2/x/(1-x))
    Gwavy_org = Gwavy

    #Find lowest point on curve
    xn=x
    dx=0.01
    dG=-1
    Gwavy2=Gwavy;
    while dG<0:
        xn=xn-dx
        if xn<0:
            xn=0
        en=F_Void_8(G,xn,Dl,Dv,ST);
        hldn=F_hld_11( en );
        Avdn=F_Avd_10( d,A,en );

        Gwavy2=(((16*Avdn**3*g*d*Dl*Dv)/(xn**2*pi**2*(1-(2*hldn-1)**2)**(1/2)))*((pi**2/(25*hldn**2))*(Wel_Frl)**(-1.023)+1))**(1/2)+50-75*exp(-(xn**2-0.97)**2/xn/(1-xn));
        dG=Gwavy2-Gwavy
        if Gwavy2<Gwavy and x>0.3:
            Gwavy=Gwavy2
    return Gwavy

def F_Gbub_26(d, A, x, e, Dl, Dv, Vl):
    g = 9.81
    #Perimeter of Interface
    Pid=F_PiD_12(e);

    #Dimensionless Cross Sectional Area Occupied by Vapor Phase
    Avd=F_Avd_10( d,A,e );

    #Dimensionless Cross Sectional Area Occupied by Liquid Phase
    Ald=F_Ald_9( d,A,e );

    #Transition Boundary - Intermittent to Bubbly Flow (I-B)
    Gbub=((256*Avd*Ald**2*d**(1.25)*Dl*(Dl-Dv)*g)/(0.3164*(1-x)**(1.75)*pi**2*Pid*Vl**(0.25)))**(1/1.75);
    return Gbub

def F_Gbub_22con(d, A, x, e, Dl, Dv, Vl):
    g = 9.81

    #Perimeter of Interface
    # Pid=F_PiD_12(e)

    #Dimensionless Cross Sectional Area Occupied by Vapor Phase
    # Avd=F_Avd_10( d,A,e )

    #Dimensionless Cross Sectional Area Occupied by Liquid Phase
    # Ald=F_Ald_9( d,A,e )

    #Transition Boundary - Intermittent to Bubbly Flow (I-B)
    #Gbub=((256*Avd*Ald**2*d**(1.25)*Dl*(Dl-Dv)*g)/(0.3164*(1-x)**(1.75)*pi**2*Pid*Vl**(0.25)))**(1/1.75);
    Gbub=F_Gbub_26(d,A,x,e,Dl,Dv,Vl);
    return Gbub

def F_Gmist_19con(G, d, A, x, e, Dl, Dv, ST):
    # G = 300
    g = 9.81
    # G=1; # G as fake input to weber and Froude, will be divided away.
    #Dimensionless Cross Sectional Area Occupied by Liquid Phase
    Ald=F_Ald_9( d,A,e );
    Avd=F_Avd_10( d,A,e );

    #Liquid Froude Number
    #Frl=F_FroudeLiquid_15( G,d,Dl );

    #Liquid Weber Number for Flow Pattern
    #Wel=F_Weber_16_21( G,d,Dl,ST);

    Wel_Frl=g*d**2*Dl/ST
    Fac=(1.138+2*log10(pi/1.5/Ald))**(-2)
    Gmist = (7680*(Avd**2)*g*d*Dl*Dv*(Wel_Frl)**(-1)/(x**2*(pi**2)*Fac))**0.5

    Gmist_org=Gmist

    # Find lowest point on curve
    xn=x
    dx=0.01
    dG=-1
    while dG.real<0:
        xn=xn-dx
        if xn<0:
            xn=0
        en=F_Void_8(G,xn,Dl,Dv,ST);
        Aldn=F_Ald_9( d,A,en );
        Avdn=F_Avd_10( d,A,en );
        Facn=(1.138+2*log10(pi/1.5/Aldn))**(-2);

        Gmist2 = (7680*(Avdn**2)*g*d*Dl*Dv*(Wel_Frl)**(-1)/(xn**2*(pi**2)*Facn))**0.5;
        dG=Gmist2-Gmist
        if Gmist2.real<Gmist.real and x>0.3:
            Gmist=Gmist2;

    #Gmist=Gmist_org;
    #--------evaporation----
    # qcrit=F_CHF_23(Dl,Dv,Hl,Hv,ST);
    # Gdry=F_Gdry_19(q,d,x,Dl,Dv,Hl,Hv,ST);
    # #Transition Boundary - Dryout to Mist Flow (D-M) (Thome et al.)
    # Gmist=((1/0.502)*(log(0.61/x)+0.57)*(d/(Dv*ST))**(-0.16)*(1/(g*d*Dv*(Dl-Dv)))**(-0.15)*(Dv/Dl)**(0.09)*(q/qcrit)**(-0.72))**(1.613);
    # if Gmist<=Gdry
    #     Gmist=Gdry;
    return Gmist

def F_PiD_12(e):
    StratAngle = F_StratAngle_13(e)
    PiD=sin((2*pi-StratAngle)/2)
    return PiD

def F_FilmAngle_8130con(G, d, A, x, e, Dl, Dv, Vl, ST):
    Gstrat=F_Gstrat_17con(d,A,x,e,Dl,Dv,Vl)
    Gwavy=F_Gwavy_16con(G,d,A,x,e,Dl,Dv,ST)
    AngleStrat=F_StratAngle_13(e)
    FilmAngle=AngleStrat*((Gwavy-G)/(Gwavy-Gstrat))**0.5
    return FilmAngle

def F_RElo_38(G, d, Vl):
    RElo=G*d/Vl
    return RElo

def F_flo_37(G, d, Vl):
    RElo=F_RElo_38(G,d,Vl);
    flo=BlasiusCorrelation(RElo)[0]/4;
    return flo

def F_DPlo_36(G, d, Dl, Vl):
    flo=F_flo_37(G,d,Vl)
    DPlo=(4*flo/d)*(G**2/(2*Dl))
    return DPlo

def F_DPbub_56(G, d, x, e, Dl, Dv, Vl, Vv, ST):
    DPa=F_DPannu_29(G,d,x,e,Dl,Dv,Vv,ST);
    xia=F_xia_18(Dl,Dv,Vl,Vv);
    eia=F_Void_8(G,xia,Dl,Dv,ST);
    DPlo=F_DPlo_36(G,d,Dl,Vl);
    DPb=DPlo*(1-e/eia)+DPa*(e/eia);
    return DPb

def F_DPannu_29(G, d, x, e, Dl, Dv, Vv, ST):
    #Mean Vapor Phase Velocity
    uv=F_uv_31(G,x,e,Dv);
    #2-Phase Friction Factor for Annular Flow
    fa= F_fa_30( G,d,x,e,Dl,Vv,ST )
    #Frictional Pressure Drop for Annular Flow
    DPa=(4*fa/d)*((Dv*uv**2)/2)
    return DPa

def F_uv_31(G, x, e, Dv):
    uv=G*x/(Dv*e)
    return uv

def F_Rev_32(G, d, x, e, Vv):
    Rev=G*x*d/(Vv*e)
    return Rev

def F_ul_34(G, x, e, Dl):
    ul=G*(1-x)/(Dl*(1-e))
    return ul

def F_WeberLiquid_33(G, d, x, e, Dl, ST):
    ul = F_ul_34(G, x, e, Dl);
    Wel = Dl*ul**2*d/ST;
    return Wel

def F_fa_30(G, d, x, e, Dl, Vv, ST):
    Rev=F_Rev_32(G,d,x,e,Vv);
    Wel=F_WeberLiquid_33(G,d,x,e,Dl,ST);
    fa=3.128*Rev**(-0.454)*Wel**(-0.0308);
    return fa

def F_DPmist_45(G, d, x, Dl, Dv, Vl, Vv):
    #Homogeneous Dynamic Viscosity (Ciccitti)
    VH=Vl*(1-x)+Vv*x;
    #Mist Flow Reynolds Number
    Rem=(G*d)/VH;
    #Mist Flow Friction Factor
    fm=91.2/Rem**(0.832);
    #Homogeneous Void Fraction
    eH=(1+((1-x)/x)*(Dv/Dl))**(-1);
    #Homogeneous Density
    DH=Dl*(1-eH)+Dv*eH;
    #Frictional Pressure Drop for Mist Flow
    DPm=(4*fm/d)*(G**2/(2*DH));
    return DPm

def F_fv_43(G, d, x, e, Vv):
    REv=F_Rev_32(G,d,x,e,Vv)
    fv=BlasiusCorrelation(REv)[0]/4
    return fv

def F_fs_53(G, d, x, e, Dl, Vv, ST):
    StratAngle = F_StratAngle_13(e)/(2*pi)
    fa=F_fa_30(G, d, x, e, Dl, Vv, ST)
    fv=F_fv_43(G,d,x,e,Vv)
    fs=StratAngle*fv+(1-StratAngle)*fa
    return fs

def F_DPstrat_52(G, d, x, e, Dl, Dv, Vl, Vv, ST):
    uv=F_uv_31(G, x, e, Dv);
    fs=F_fs_53(G, d, x, e, Dl, Vv, ST);
    DPs=4*fs*Dv*uv**2/(2*d);
    xia=F_xia_18(Dl,Dv,Vl,Vv);
    eia=F_Void_8(G,xia,Dl,Dv,ST);
    if x < xia:
        DPlo=F_DPlo_36(G,d,Dl,Vl);
        DPs=DPlo*(1-e/eia)+DPs*(e/eia);
    return DPs

def F_DPsw_39(G, d, A, x, e, Dl, Dv, Vl, Vv, ST):
    #Mean Vapor Phase Velocity
    uv=F_uv_31(G,x,e,Dv);

    fsw=F_fsw_40( G,d,A,x,e,Dl,Dv,Vl,Vv,ST );

    #Frictional Pressure Drop for Stratified-Wavy Flow
    #DPsw=4*fsw*((l*Dv*Vv**2)/(2*d)); bug in Joao's code
    DPsw=4*fsw*((Dv*uv**2)/(2*d));
    return DPsw

def F_DryAngle_42( G,d,A,x,e,Dl,Dv,Vl,Vv,ST ):
    StratAngle=F_StratAngle_13(e);
    Gwavy=F_Gwavy_14(G,d,A,x,e,Dl,Dv,ST);
    Gstrat= F_Gstrat_17(G,d,A,x,e,Dl,Dv,Vl,Vv,ST) ;
    DryAngle=StratAngle*((Gwavy-G)/(Gwavy-Gstrat))**0.61;
    return DryAngle



def F_fsw_40(G, d, A, x, e, Dl, Dv, Vl, Vv, ST):
    DryAngle=F_DryAngle_42( G,d,A,x,e,Dl,Dv,Vl,Vv,ST )/(2*pi);
    fa=F_fa_30(G,d,x,e,Dl,Vv,ST);
    fv=F_fv_43(G,d,x,e,Vv);
    fsw=DryAngle**0.02*fv+(1-DryAngle)**0.02*fa;
    return fsw

def F_DPslug_int_35(G, d, x, e, Dl, Dv, Vl, Vv, ST):
    #Delta P for slug and intermittent flow
    DPlo=F_DPlo_36(G,d,Dl,Vl);
    DPa=F_DPannu_29(G,d,x,e,Dl,Dv,Vv,ST);
    xia=F_xia_18(Dl,Dv,Vl,Vv);
    eia=F_Void_8(G,xia,Dl,Dv,ST);

    DPsi=DPlo*(1-(e/eia))+DPa*(e/eia);
    return DPsi

def F_DPmomentum_28(G, q, A, Ph, x, e, Dl, Dv, Hl, Hv, ST):
    l=0.01;
    Q=q*Ph*l;
    MF=G*A;
    xin=x;
    ein=e;
    xout=Q/(MF*(Hv-Hl))+xin;
    eout=F_Void_8(G,xout,Dl,Dv,ST);
    dPm=G**2*((((1-xin)**2/Dl*(1-ein))+(xin**2/Dv*ein))-(((1-xout)**2/Dl*(1-eout))+(xout**2/Dv*eout)))/l;
    return dPm

def F_HTC_8123con(G, q, d, A, x, e, FilmAngle, Dl, Dv, CPl, Kl, Hl, Hv, Vl, ST):
    #Fluid,P,H,G,q,d,A,Ph,Ep,Angle,x,e,Angle_dry,Dl,Dv,Cpl,Kl,Kv,Hl,Hv,Vl,Vv,ST
    #htc
    htcAxial=F_htcAxial_8132con(G,d,A,x,e,Dl,Dv,CPl,Kl,Vl,ST);
    #htf
    htcFilm=F_htcFilm_8143con(q,d,Dl,Dv,Kl,Vl,Hl,Hv);
    #htp
    htc=(htcFilm*FilmAngle+(2*pi-FilmAngle)*htcAxial)/(2*pi);
    return htc

def F_htcFilm_8143con(q, d, Dl, Dv, Kl, Vl, Hl, Hv):
    #Fluid,P,H,G,q,d,A,Ph,Ep,Angle,x,e,Angle_dry,Dl,Dv,Kl,Kv,Hl,Hv,Vl,Vv,ST
    g=9.81;
    htf=0.655*(Dl*(Dl-Dv)*g*(Hv-Hl)*Kl**3/(Vl*d*q))**(1/3);
    return htf

def F_htcAxial_8132con(G, d, A, x, e, Dl, Dv, CPl, Kl, Vl, ST):
    c=0.003
    n=0.74
    m=0.5
    dfilm=F_dfilm_8135con(G,d,A,x,e,Dl,Dv,Vl,ST);
    Rel=F_ReL_8133con(G,x,e,Vl,dfilm);
    Pr=F_Prandtx( CPl,Kl,Vl) ;
    Gstrat=F_Gstrat_17con(d,A,x,e,Dl,Dv,Vl);
    fi=F_fi_8140con(G,d,A,x,e,Dl,Dv,Vl,ST);
    if G<Gstrat:
        fi=fi*G/Gstrat
    htc=c*Rel**n*Pr**m*Kl/dfilm*fi;
    return htc

def F_fi_8140con(G, d, A, x, e, Dl, Dv, Vl, ST):
    g=9.81;
    dfilm=F_dfilm_8135con(G,d,A,x,e,Dl,Dv,Vl,ST);
    uv=G*x/(Dv*e);
    ul=G*(1-x)/(Dl*(1-e));
    fi=1+(uv/ul)**0.5*((Dl-Dv)*g*dfilm**2/ST)**0.25;
    return fi

def F_dfilm_8135con(G, d, A, x, e, Dl, Dv, Vl, ST):
    Al=A*(1-e)
    FilmAngle=F_FilmAngle_8130con(G,d,A,x,e,Dl,Dv,Vl,ST);
    dfilm=(d-(d**2-Al*8/(2*pi-FilmAngle))**0.5)/2;
    # dAl=1;
    # dfilm=0;
    # while dAl>1e-20;
    # Al_guess=((2*pi-FilmAngle)/8)*(d**2-(d-2*dfilm)**2);
    # dAl=Al-Al_guess;
    # dfilm=dfilm+1e-8;

    dfilm=(d/2)-sqrt((d/2)**2-(2*Al)/(2*pi-FilmAngle));
    if dfilm>d/2:
        dfilm=d/2
    return dfilm

def F_ReL_8133con(G, x, e, Vl, df):
    Rel=4*G*(1-x)*df/(1-e)/Vl
    return Rel

def F_Prandtx(CPx, Kx, Vx):
    #Liquid or vapor Prandtl number
    PR=(CPx*Vx)/Kx;
    return PR

def F_Gwavy_14(G, d, A, x, e, Dl, Dv, ST):
    #Kattan�Thome�Favrat criterion
    g=9.81;

    #Dimensionless Cross Sectional Area Occupied by Vapor Phase
    Avd=F_Avd_10( d,A,e );

    #Dimensionless Vertical Height of the Liquid
    hld=F_hld_11( e );

    #Liquid Froude Number
    #Frl=F_FroudeLiquid_15( G,d,Dl );

    #Liquid Weber Number for Flow Pattern

    #Wel=F_Weber_16_21( G,d,Dl,ST);

    Wel_Frl=g*d**2*Dl/ST; # (We/Fr)l

    #Transition Boundary - Stratified-Wavy to Intermittent and Annular (SW-I/A)
    Gwavy=(((16*Avd**3*g*d*Dl*Dv)/(x**2*pi**2*(1-(2*hld-1)**2)**(1/2)))*(((pi**2)/(25*hld**2))*(Wel_Frl)**-1+1))**(1/2)+50;
    return Gwavy

def F_Gstrat_17(G, d, A, x, e, Dl, Dv, Vl, Vv, ST):
    #Fluid,P,H,G,q,d,A,Ph,Ep,Angle,x,e,Angle_dry,Dl,Dv,Hl,Hv,Vl,Vv
    g=9.81;
    xia=F_xia_18(Dl,Dv,Vl,Vv);
    eia=F_Void_8(G,xia,Dl,Dv,ST);
    #Dimensionless Cross Sectional Area Occupied by Vapor Phase
    Avd=F_Avd_10( d,A,e );
    Avd_xia=F_Avd_10( d,A,eia );
    #Dimensionless Cross Sectional Area Occupied by Liquid Phase
    Ald=F_Ald_9( d,A,e );
    Ald_xia=F_Ald_9( d,A,eia );

    #Transition Boundary - Stratified to Stratified-Wavy (S-SW)
    if x<xia:
        Gstrat=((226.3**2*Ald_xia*Avd_xia**2*Dv*(Dl-Dv)*Vl*g)/(xia**2*(1-xia)*pi**3))**(1/3);
    else:
        Gstrat=((226.3**2*Ald*Avd**2*Dv*(Dl-Dv)*Vl*g)/(x**2*(1-x)*pi**3))**(1/3);
    return Gstrat

def F_Gmist_24(q,d,x,Dl,Dv,Hl,Hv,ST):
    g=9.81
    qcrit=F_CHF_23(Dl,Dv,Hl,Hv,ST);
    Gdry=F_Gdry_19(q,d,x,Dl,Dv,Hl,Hv,ST);
    #Transition Boundary - Dryout to Mist Flow (D-M) (Thome et al.)
    Gmist=((1/0.502)*(log(0.61/x)+0.57)*(d/(Dv*ST))**(-0.16)*(1/(g*d*Dv*(Dl-Dv)))**(-0.15)*(Dv/Dl)**(0.09)*(q/qcrit)**(-0.72))**(1.613);
    if Gmist<=Gdry:
        Gmist=Gdry
    return Gmist

def F_FroudeVapor_22( G,d,Dl,Dv ):
    # Froude vapor by Mori et al.
    g=9.81;
    Frv=G**2/(Dv*(Dl-Dv)*g*d);
    return Frv

def F_FroudeLiquid_15( G,d,Dl ):
    g=9.81;
    Frl=G**2/(Dl**2*g*d);
    return Frl

def F_CHF_23(Dl, Dv, Hl, Hv, ST):
    g=9.81
    #Critical Heat Flux (By Kutateladze)
    qcrit=0.131*Dv**(0.5)*(Hv-Hl)*(g*ST*(Dl-Dv))**(0.25);
    return qcrit

def F_Gdry_19(q, d, x, Dl, Dv, Hl, Hv, ST):
    g=9.81;
    qcrit=F_CHF_23(Dl,Dv,Hl,Hv,ST)
    Gdry=((1/0.236)*(log(0.58/x)+0.52)*(d/(Dv*ST))**(-0.17)*(1/(g*d*Dv*(Dl-Dv)))**(-0.17)*(Dv/Dl)**(-0.25)*(q/qcrit)**(-0.27))**(1.471);
    return Gdry

def F_Weber_16_21(G, d, Dx, ST):
    #Weber nummer. Dx is either the liquid or the vapor density
    We=(G**2*d)/(Dx*ST);
    return We

def F_xde_25(G, q, d, Dl, Dv, Hl, Hv, ST):
    #Dry out completion vapor quality
    #Fluid,P,H,G,q,d,A,Ph,Ep,Angle,x,e,Angle_dry,Dl,Dv,Hl,Hv,Vl,Vv,ST
    qcrit=F_CHF_23(Dl,Dv,Hl,Hv,ST);
    Wev=F_Weber_16_21( G,d,Dv,ST );
    Frv=F_FroudeVapor_22( G,d,Dl,Dv );
    xde=0.61*exp(0.57-0.502*Wev**(0.16)*Frv**(0.15)*(Dv/Dl)**(-0.09)*(q/qcrit)**(0.72));
    return xde

def F_xdi_20(G, q, d, Dl, Dv, Hl, Hv, ST):
    #Dry out inception vapor quality
    qcrit = F_CHF_23(Dl, Dv, Hl, Hv, ST)
    Wev = F_Weber_16_21(G, d, Dv, ST )
    Frv = F_FroudeVapor_22(G, d, Dl, Dv)
    xdi=0.58*exp(0.52-0.236*Wev**(0.17)*Frv**(0.17)**(Dv/Dl)**(0.25)*(q/qcrit)**(0.27))
    return xdi

def F_DPslug_sw_44(G,d,A,x,e,Dl,Dv,Vl,Vv,ST):
    # Delta P for slug and stritified wavy flow
    DPlo=F_DPlo_36(G,d,Dl,Vl);
    DPsw=F_DPsw_39(G,d,A,x,e,Dl,Dv,Vl,Vv,ST);
    xia=F_xia_18(Dl,Dv,Vl,Vv);
    eia=F_Void_8(G,xia,Dl,Dv,ST);

    DPssw=DPlo*(1-(e/eia))+DPsw*(e/eia);
    return DPssw

def F_DPdry_51(G, q, d, A, x, Dl, Dv, Hl, Hv, Vl, Vv, ST):
    xia = F_xia_18(Dl,Dv,Vl,Vv);
    xdi = F_xdi_20( G,q,d,Dl,Dv,Hl,Hv,ST);
    xde = F_xde_25( G,q,d,Dl,Dv,Hl,Hv,ST);
    edi = F_Void_8(G,xdi,Dl,Dv,ST);
    Gwavy = F_Gwavy_14(G,d,A,xdi,edi,Dl,Dv,ST);
    DPm = F_DPmist_45(G,d,xde,Dl,Dv,Vl,Vv);
    if G >= Gwavy and xdi < xia:
        DP=F_DPbub_56(G,d,xdi,edi,Dl,Dv,Vl,Vv,ST ) ;
    elif G >= Gwavy and xdi >= xia:
        DP=F_DPannu_29(G,d,xdi,edi,Dl,Dv,Vv,ST);
    else:
        DP=F_DPsw_39(G,d,A,xdi,edi,Dl,Dv,Vl,Vv,ST);
    DPd=DP-(((x-xdi)/(xde-xdi))*(DP-DPm));
    return DPd

def F_PiD_12(e):
    StratAngle=F_StratAngle_13(e);
    PiD=sin((2*pi-StratAngle)/2);
    return PiD

def F_HTCannu_int_bub_slug(p, G, q, d, A, x, e, CPl, CPv, Dl, Dv, Kl, Kv, Vl, Vv, ST):
    DryAngle=0;
    htp=F_HTCtp_1( p,G,q,d,A,x,e,CPl,CPv,Dl,Dv,Kl,Kv,Vl,Vv,ST,DryAngle);
    return htp

def F_HTCstrat( p,G,q,d,A,x,e,CPl,CPv,Dl,Dv,Kl,Kv,Vl,Vv,ST):
    DryAngle=F_StratAngle_13(e);
    htp=F_HTCtp_1( p,G,q,d,A,x,e,CPl,CPv,Dl,Dv,Kl,Kv,Vl,Vv,ST,DryAngle);
    return htp

def F_HTCtp_1045evap( fluid,G,q,d,A,x,e,CPl,Cpv,Kl,Kv,Vl,Vv,DryAngle):
    # hvap=F_hvap_1049evap(G,d,x,e,Cpv,Kv,Vv);
    # hwet=F_hwet_1046evap( fluid,G,q,d,A,x,e,CPl,Kl,Vl,DryAngle);
    hvap = hwet = 0.2
    htp=(DryAngle*hvap+(2*pi-DryAngle)*hwet)/(2*pi);
    return htp

def F_HTCsw( p,G,q,d,A,x,e,CPl,CPv,Dl,Dv,Kl,Kv,Vl,Vv,ST):
    DryAngle=F_DryAngle_42( G,d,A,x,e,Dl,Dv,Vl,Vv,ST );

    htp=F_HTCtp_1( p,G,q,d,A,x,e,CPl,CPv,Dl,Dv,Kl,Kv,Vl,Vv,ST,DryAngle);
    return htp

def F_HTCslug_sw( p,G,q,d,A,x,e,CPl,CPv,Dl,Dv,Kl,Kv,Vl,Vv,ST):
    xia=F_xia_18(Dl,Dv,Vl,Vv);
    DryAngle=F_DryAngle_42( G,d,A,x,e,Dl,Dv,Vl,Vv,ST )*x/xia;
    htp=F_HTCtp_1( p,G,q,d,A,x,e,CPl,CPv,Dl,Dv,Kl,Kv,Vl,Vv,ST,DryAngle);
    return htp

def F_HTCmist_14( G,d,x,CPv,Dl,Dv,Kv,Vv ):
    Reh=(G*d/Vv)*(x+(Dv/Dl)*(1-x));   #equation 15
    Y=1-0.1*(((Dl/Dv)-1)*(1-x))**0.4;  #equation 16
    PRv=F_Prandtl( CPv,Kv,Vv);
    HTCmist=2e-8*Reh**1.97*PRv**(1.06)*Y**(-1.83)*Kv/d;
    return HTCmist

def F_HTCdry_17( p,G,q,d,A,x,CPl,CPv,Dl,Dv,Hl,Hv,Kl,Kv,Vl,Vv,ST):
    xde=F_xde_25( G,q,d,Dl,Dv,Hl,Hv,ST);
    xdi=F_xdi_20( G,q,d,Dl,Dv,Hl,Hv,ST);
    edi=F_Void_8(G,xdi,Dl,Dv,ST);
    HTCxde=F_HTCmist_14( G,d,xde,CPv,Dl,Dv,Kv,Vv );
    Gwavy=F_Gwavy_14(G,d,A,xdi,edi,Dl,Dv,ST);
    if G>=Gwavy:
        DryAngle=0;
    else:
        DryAngle=F_DryAngle_42( G,d,A,xdi,edi,Dl,Dv,Vl,Vv,ST );

    HTCxdi=F_HTCtp_1( p,G,q,d,A,xdi,edi,CPl,CPv,Dl,Dv,Kl,Kv,Vl,Vv,ST,DryAngle);
    HTCdry=HTCxdi-((x-xdi)/(xde-xdi))*(HTCxdi-HTCxde);
    return HTCdry



def F_Prandtl(CPx, Kx, Vx):
    #Liquid or vapor Prandtl number
    PRx=(CPx*Vx)/Kx;
    return PRx

def F_hnb_8( p,q):
    #CO2 Critical Pressure
    pcrit=73.773;
    #CO2 Molar Mass
    M=44.0098;
    #Reduced Pressure
    pr=p/pcrit;
    #Nucleate Boiling HTC
    hnb=131*pr**(-0.0063)*(-log10(pr))**(-0.55)*M**(-0.5)*q**(0.58);
    return hnb

def F_delta_11(d,A,e,DryAngle):
    delta=(d/2)-((d/2)^2-((2*A*(1-e))/(2*pi-DryAngle)))**0.5;
    if delta>d/2:
        delta=d/2;
    return delta

def F_REd_13(G,d,A,x,e,Vl,DryAngle ):
    delta= F_delta_11(d,A,e,DryAngle);
    Red=4*G*(1-x)*delta/(Vl*(1-e));
    return Red

def F_hcb_12(G, d, A, x, e, CPl, Kl, Vl, DryAngle):
    Red=F_REd_13(G,d,A,x,e,Vl,DryAngle );
    PRl=F_Prandtl( CPl,Kl,Vl);
    delta=F_delta_11(d,A,e,DryAngle);
    hcb=0.0133*Red**0.69*PRl**0.4*Kl/delta;
    return hcb

def F_hv_5(G, d, x, e, CPv, Kv, Vv):
    Rev=F_Rev_32(G,d,x,e,Vv);
    Prv=F_Prandtl(CPv,Kv,Vv);
    hv=0.023*Rev**0.8*Prv**0.4*Kv/d;
    return hv

def F_Supp_9_10( G,d,A,x,e,Dl,Dv,Vl,Vv,ST,DryAngle):
    #Boiling Suppression Factor
    xia=F_xia_18(Dl,Dv,Vl,Vv);
    eia=F_Void_8(G,xia,Dl,Dv,ST);
    delta=F_delta_11(d,A,e,DryAngle);
    delta_ia=F_delta_11(d,A,eia,DryAngle);
    if x<xia:
        S=1;
    else:
        S=1-1.14*(min(d, 0.00753)/0.00753)**2*(1-(delta/delta_ia))**2.2;
    return S

def F_hwet_7(p, G, q, d, A, x, e, CPl, Dl, Dv, Kl, Vl, Vv, ST, DryAngle):
    S=F_Supp_9_10( G,d,A,x,e,Dl,Dv,Vl,Vv,ST,DryAngle );
    hnb=F_hnb_8( p,q);
    hcb=F_hcb_12( G,d,A,x,e,CPl,Kl,Vl,DryAngle );
    hwet=((S*hnb)**3+hcb**3)**(1/3);
    return hwet

def F_HTCtp_1(p, G, q, d, A, x, e, CPl, CPv, Dl, Dv, Kl, Kv, Vl, Vv, ST, DryAngle):
    hv=F_hv_5( G,d,x,e,CPv,Kv,Vv );
    hwet=F_hwet_7( p,G,q,d,A,x,e,CPl,Dl,Dv,Kl,Vl,Vv,ST,DryAngle  );
    htp=(DryAngle*hv+(2*pi-DryAngle)*hwet)/(2*pi);
    return htp
