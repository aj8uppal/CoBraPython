# Pressure drop and HTC function.  Determines which pressure drop and heattransfg=er models to use
# based on the PH situation.
# Inputs: Fluid, P (bar), H (J/kg), MFLX (kg/m2s), Dh (m), Ep (m) & Angle (Degree, Super heating ('C))
# Output: dP (bar/m), HTC (W/m2K), Vapor quality VQ (-), relative mass
# (kg/m), State, Temperature ('C)
import sys
sys.path.append('..')
from FluidDynamicEquations import *
sys.path.append('../FluidDynamicEquations')
from Correlations import *
# sys.path.append('../../REFPROP')
# from refprop import RefPropInterface
from math import sin, pi
from time import time

sind = lambda x: sin(x)*180/pi

def dPandHTC(Fluid,P,H,MFLX,HFLX,Dh,A,Ph,Ep,Angle,SH, refpropm):
    # RPI = RefPropInterface()
    # refpropm = RPI.refpropm
    dP=-1;
    Fd=-1;
    HTC=-1;
    VQ=-1;
    Void=-1;
    State='null';
    g=9.81;
    rm=-1;

    # ver=version;
    # ver=str2num(ver(end-5:end-2));
    # if ver>2011;
    times = {1: 0, 2: 0}
    # CritP=refpropm('P','C',0,' ',0,Fluid)*1e-2;
    #
    # if Fluid[-3:].lower() == 'mix':
    #     TripP = 0
    # else:
    #     TripP=refpropm('P','R',0,' ',0,Fluid)*1e-2;
    TripP = 5.1796
    CritP = 73.7730


    #if fluid is co2 CritP=73.7; TripP=5.2;
    #else CritP=1e5;TripP=0;

    #Define fluid state.

    if P >= CritP:
        State='supercritical'
    elif P <= TripP:
        Hsatgas=refpropm('H','P',P*1e2,'Q',1,Fluid)
        if H < Hsatgas:
            State='sublimation'
        else:
            State='gas'
    else:
        Hsatliq=refpropm('H','P',P*1e2,'Q',0,Fluid);
        Hsatgas=refpropm('H','P',P*1e2,'Q',1,Fluid);
        Cp=refpropm('C','P',P*1e2,'Q',0,Fluid);
        Tsh=(H - Hsatliq)/Cp;
        if H <= Hsatliq:
            State='liquid'
        elif H>=Hsatgas:
            State='gas'
        elif Tsh<=SH:
            State='superheated'
        else:
            State='2phase'
            VQ=(H-Hsatliq)/(Hsatgas-Hsatliq)

    T=refpropm('T','P',P*1e2,'H',H,Fluid)-273.15;
    # times[1]+=time()-start
    # start = time()
    if State == 'supercritical':

        dP, rm, none = DarcyWeisbach(Fluid,P,H,MFLX,Dh,A,Ep, refpropm)

        HTC = DittusBoelter( Fluid,P,H,MFLX,HFLX,Dh, refpropm)
        VQ=0
    elif State == 'liquid':
        dP, rm, State = DarcyWeisbach(Fluid,P,H,MFLX,Dh,A,Ep, refpropm);
        HTC = DittusBoelter( Fluid,P,H,MFLX,HFLX,Dh, refpropm);
        VQ=0;
    elif State == 'superheated':
        dP, rm, none = DarcyWeisbach(Fluid,P,Hsatliq-1,MFLX,Dh,A,Ep, refpropm);
        HTC = DittusBoelter( Fluid,P,Hsatliq-1,MFLX,HFLX,Dh, refpropm);
        VQ=0;
        T=T+Tsh;
    elif State == 'gas':
        dP, rm, none = DarcyWeisbach(Fluid,P,H,MFLX,Dh,A,Ep, refpropm);
        HTC = DittusBoelter( Fluid,P,H,MFLX,HFLX,Dh, refpropm);
        VQ=1
    elif State == 'sublimation':
        dP=0;
        HTC=0;
        VQ=0;
        rm=0;
    elif State == '2phase':
        if HFLX>=0:
            if Angle==0:
                if Fluid == 'CO2':
                    # print("{}, {}, {}, {}, {}, {}, {}, {}".format(Fluid,P,H,MFLX,HFLX,Dh,A,Ph))
                    dP,HTC,VQ,rm,State=ThomeCorrelation_EvapHor_CO2(Fluid,P,H,MFLX,HFLX,Dh,A,Ph, refpropm);
                else:
                     #HTC=KandlikarCorrelation(Fluid,P,VQ,MFLX,HFLX,Dh);
                     #[dP,HTC,VQ,rm,State]=ThomeCorrelation_EvapHor(Fluid,P,H,MFLX,HFLX,Dh,A,Ph);
                     dP,rm=FriedelCorrelation(Fluid,P,VQ,MFLX,Dh,A, refpropm);
                     HTC=KandlikarCorrelation(Fluid,P,VQ,MFLX,HFLX,Dh, refpropm)
            elif abs(Angle)==90:
                #vertical formulas
                dP, rm = FriedelCorrelation(Fluid,P,VQ,MFLX,Dh,A, refpropm);
                HTC=KandlikarCorrelation(Fluid,P,VQ,MFLX,HFLX,Dh, refpropm);
            else:
                #inclined formulas
                dP, rm = FriedelCorrelation(Fluid,P,VQ,MFLX,Dh,A, refpropm);
                HTC=KandlikarCorrelation(Fluid,P,VQ,MFLX,HFLX,Dh, refpropm);
        elif HFLX<0:
            dP,HTC,VQ,rm,State = ThomeCorrelation_Con(Fluid,P,H,MFLX,HFLX,Dh,A,Ph, refpropm);

    # times[2]+=time()-start
    dPstat=1e-5*g*rm*sind(Angle)/A;
    dP=dP+dPstat;
    # print(times)
    return dP,HTC,VQ,rm,State,T
