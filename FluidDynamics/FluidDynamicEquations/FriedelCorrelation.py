from Correlations import F_Void_8
import sys
sys.path.append('../../REFPROP')
from refprop import RefPropInterface

#VERIFIED

def FriedelCorrelation(Fluid, P, VQ, MFLX, Dh, A):
    RPI = RefPropInterface()
    refpropm = RPI.refpropm
    #Fluid properties: (arbitrary values until refprop is fixed)
    Dliq=refpropm('D','P',P*1e2,'Q',0,Fluid);       #kg/m3
    Dvap=refpropm('D','P',P*1e2,'Q',1,Fluid);      #kg/m3
    Dtp=refpropm('D','P',P*1e2,'Q',VQ,Fluid);       #kg/m3
    Vliq=refpropm('V','P',P*1e2,'Q',0,Fluid);    #Pa*s
    Vvap=refpropm('V','P',P*1e2,'Q',1,Fluid);     #Pa.s
    Isft=refpropm('I','P',P*1e2,'Q',0,Fluid);    #N/m

    # Dliq = Dvap = Dtp = Vliq = Vvap = Isft = 2

    #Dimensionless numbers
    REliq=MFLX*(1-VQ)*Dh/Vliq    #Reynolds number of the liquid phase
    REvap=MFLX*VQ*Dh/Vvap        #Reynolds number of the vapor phase
    REliqO=MFLX*Dh/Vliq            #Reynolds number assuming liquid only
    REvapO=MFLX*Dh/Vvap           #Reynolds number assuming vapor only
    FR=MFLX**2/(9.81*Dh*Dtp**2)      #Froude number
    WE=MFLX**2*Dh/(Dtp*Isft)        #Weber number

    #Calculate friction factors  (Blasius equation used here)
    Fliq=0.079/REliq**(1/4)         #Liquid friction factor
    Fvap=0.079/REvap**(1/4)         #Vapor friction factor
    FliqO=0.079/REliqO**(1/4)         #Liquid only friction factor
    FvapO=0.079/REvapO**(1/4)         #Vapor only friction factor

    #Calculate equivilent pressure gradients
    dPliq=2*Fliq*MFLX**2*(1-VQ)**2/(Dh*Dliq)/1e5    #bar/m
    dPvap=2*Fvap*MFLX**2*VQ**2/(Dh*Dvap)/1e5    #bar/m
    dPliqO=2*FliqO*MFLX**2/(Dh*Dliq)/1e5    #bar/m
    dPvapO=2*FvapO*MFLX**2/(Dh*Dvap)/1e5    #bar/m

    E=(1-VQ)**2+VQ**2*(Dliq*FvapO/(Dvap*FliqO))
    F=VQ**0.78*(1-VQ)**0.24
    H=(Dliq/Dvap)**0.91*(Vvap/Vliq)**0.19*(1-(Vvap/Vliq))**0.7
    TPmult2=E+(3.24*F*H)/(FR**0.045*WE**0.035)   #2-Phase multiplier according to Friedel

    dP=TPmult2*dPliqO #2-phase pressure gradient bar/m

    e=F_Void_8(MFLX,VQ,Dliq,Dvap,Isft) #Void fraction
    rm=(Dvap*e+Dliq*(1-e))*A #Relative mass kg/m

    return dP, rm
