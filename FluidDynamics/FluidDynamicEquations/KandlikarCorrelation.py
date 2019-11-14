#VERIFIED
import sys
sys.path.append('../../REFPROP')
from refprop import RefPropInterface

def KandlikarCorrelation(Fluid, P, VQ, MFLX, HFLX, Dh):
    #Commented until refprop is fixed, arbitrary value for now
    #Fluid properties
    RPI = RefPropInterface()
    refpropm = RPI.refpropm
    Dliq=refpropm('D','P',P*1e2,'Q',0,Fluid)       #kg/m3
    Dvap=refpropm('D','P',P*1e2,'Q',1,Fluid)      #kg/m3
    Dtp=refpropm('D','P',P*1e2,'Q',VQ,Fluid)       #kg/m3
    Vliq=refpropm('V','P',P*1e2,'Q',0,Fluid)    #Pa*s
    Vvap=refpropm('V','P',P*1e2,'Q',1,Fluid)     #Pa.s
    Isft=refpropm('I','P',P*1e2,'Q',0,Fluid)    #N/m


    Hliq=refpropm('H','P',P*1e2,'Q',0,Fluid)       #J/kg
    Hvap=refpropm('H','P',P*1e2,'Q',1,Fluid)       #J/kg
    LAMliq=refpropm('L','P',P*1e2,'Q',0,Fluid)      # W/mK
    CPliq=refpropm('C','P',P*1e2,'Q',0,Fluid)        #Specific heat
    Dliq = Dtp = Vliq  = Isft = Hliq = LAMliq = CPliq = 2
    Dvap = Vvap = Hvap = 3
    PRliq=Vliq*CPliq/LAMliq #Prandl

    #Dimensionless numbers
    RElo=MFLX*Dh/Vliq #Reynolds number assuming liquid only

    Co=(((1-VQ)/VQ)**0.8)*((Dvap/Dliq)**0.5) #Convection number
    Bo=HFLX/(MFLX*(Hvap-Hliq)) #Boiling Number
    FRliq=MFLX**2/(9.81*Dh*Dliq**2) #Liquid Froude number
    C1con=1.136;
    C2con=-0.9;
    C3con=667.2;
    C4con=0.7;
    C5con=0

    C1nuc=0.6683;
    C2nuc=-0.2;
    C3nuc=1058;
    C4nuc=0.7;
    C5nuc=0

    if FRliq<0.04:
        C5con=0.3
        C5nuc=0.3 #C5 choice only for horizontal tubes, vertical tubes always are 0

    Ffl=1;
    NUlo=0.023*RElo**0.8*PRliq**0.4;
    ALFAlo=NUlo*LAMliq/Dh;
    ALFAcon=ALFAlo*(C1con*Co**C2con*(25*FRliq)**C5con+C3con*Bo**C4con*Ffl);
    ALFAnuc=ALFAlo*(C1nuc*Co**C2nuc*(25*FRliq)**C5nuc+C3nuc*Bo**C4nuc*Ffl);
    ALFAtp=ALFAnuc;
    return ALFAtp
