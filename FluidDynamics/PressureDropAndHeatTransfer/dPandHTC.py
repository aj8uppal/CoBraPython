# Pressure drop and HTC function.  Determines which pressure drop and heattransfg=er models to use
# based on the PH situation.
# Inputs: Fluid, P (bar), H (J/kg), MFLX (kg/m2s), Dh (m), Ep (m) & Angle (Degree, Super heating ('C))
# Output: dP (bar/m), HTC (W/m2K), Vapor quality VQ (-), relative mass
# (kg/m), State, Temperature ('C)
#

import sys
sys.path.append('../FluidDynamicEquations')
from Correlations import *
sys.path.append('../../REFPROP')
from refprop import RefPropInterface

def dPandHTC(Fluid,P,H,MFLX,HFLX,Dh,A,Ph,Ep,Angle,SH):
    RPI = RefPropInterface()
    refpropm = RPI.refpropm
    dP=-1;
    Fd=-1;
    HTC=-1;
    VQ=-1;
    Void=-1;
    State='null';
    g=9.81;
    rm=-1;

    ver=version;
    ver=str2num(ver(end-5:end-2));
    if ver>2011;
        CritP=refpropm('P','C',0,' ',0,Fluid)*1e-2;
        if strcmp(Fluid(end-2:end),'MIX')==1 || strcmp(Fluid(end-2:end),'Mix')==1 ||strcmp(Fluid(end-2:end),'mix')==1;
            TripP=0;
        else
            TripP=refpropm('P','R',0,' ',0,Fluid)*1e-2;
        end
    elseif strcmp(Fluid,'CO2')==1;
        CritP=73.7; TripP=5.2;
    else CritP=1e5;TripP=0;
    end




    % Define fluid state.

    if P>=CritP;
        State='supercritical';
    elseif P<=TripP;
        Hsatgas=refpropm('H','P',P*1e2,'Q',1,Fluid);
        if H<Hsatgas;
        State='sublimation';
        else
            State='gas';
        end
    else
        Hsatliq=refpropm('H','P',P*1e2,'Q',0,Fluid);
        Hsatgas=refpropm('H','P',P*1e2,'Q',1,Fluid);
        Cp=refpropm('C','P',P*1e2,'Q',0,Fluid);
        Tsh=(H-Hsatliq)/Cp;

        if H<=Hsatliq;
            State='liquid';
        elseif H>=Hsatgas;
            State='gas';
        elseif Tsh<=SH
            State='superheated';
        else State='2phase';
            VQ=(H-Hsatliq)/(Hsatgas-Hsatliq);
        end
    end
    % if strcmp(Fluid,'ASML')==1
    %     State='2phase';
    % end

    T=refpropm('T','P',P*1e2,'H',H,Fluid)-273.15;

    if strcmp(State,'supercritical')==1;
        [dP,rm,none]=DarcyWeisbach(Fluid,P,H,MFLX,Dh,A,Ep);
        HTC=DittusBoelter( Fluid,P,H,MFLX,HFLX,Dh);
        VQ=0;
    elseif strcmp(State,'liquid')==1;
        [dP,rm,State]=DarcyWeisbach(Fluid,P,H,MFLX,Dh,A,Ep);
        HTC=DittusBoelter( Fluid,P,H,MFLX,HFLX,Dh);
        VQ=0;
    elseif strcmp(State,'superheated')==1;
        [dP,rm,none]=DarcyWeisbach(Fluid,P,Hsatliq-1,MFLX,Dh,A,Ep);
        HTC=DittusBoelter( Fluid,P,Hsatliq-1,MFLX,HFLX,Dh);
        VQ=0;
        T=T+Tsh;
    elseif strcmp(State,'gas')==1;
        [dP,rm,none]=DarcyWeisbach(Fluid,P,H,MFLX,Dh,A,Ep);
        HTC=DittusBoelter( Fluid,P,H,MFLX,HFLX,Dh);
        VQ=1;
    elseif strcmp(State,'sublimation')==1;
        dP=0;
        HTC=0;
        VQ=0;
        rm=0;
    elseif strcmp(State,'2phase')==1;

            if HFLX>=0;
                 if Angle==0;
                    if strcmp(Fluid,'CO2')==1;
                     [dP,HTC,VQ,rm,State]=ThomeCorrelation_EvapHor_CO2(Fluid,P,H,MFLX,HFLX,Dh,A,Ph);
                    else
                     %HTC=KandlikarCorrelation(Fluid,P,VQ,MFLX,HFLX,Dh);
                     %[dP,HTC,VQ,rm,State]=ThomeCorrelation_EvapHor(Fluid,P,H,MFLX,HFLX,Dh,A,Ph);
                     [dP,rm]=FriedelCorrelation(Fluid,P,VQ,MFLX,Dh,A);
                     HTC=KandlikarCorrelation(Fluid,P,VQ,MFLX,HFLX,Dh);
                    end
                elseif Angle==90 || Angle==-90
                % vertical formulas
                    [dP,rm]=FriedelCorrelation(Fluid,P,VQ,MFLX,Dh,A);
                    HTC=KandlikarCorrelation(Fluid,P,VQ,MFLX,HFLX,Dh);
                else
                % inclined formulas
                     [dP,rm]=FriedelCorrelation(Fluid,P,VQ,MFLX,Dh,A);
                     HTC=KandlikarCorrelation(Fluid,P,VQ,MFLX,HFLX,Dh);
                end
            elseif HFLX<0;
                [dP,HTC,VQ,rm,State]=ThomeCorrelation_Con(Fluid,P,H,MFLX,HFLX,Dh,A,Ph);
            end

     end


    dPstat=1e-5*g*rm*sind(Angle)/A;
     dP=dP+dPstat;

    %  dP=real(dP);
    %  HTC=real(HTC);
    %  VQ=real(VQ);
    %  rm=real(rm);
    %  T=real(T);

    return dP,HTC,VQ,rm,State,T
