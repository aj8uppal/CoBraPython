from .Correlations import F_Void_8

#VERIFIED

def ChisholmCorrelation(Fluid, P, VQ, MFLX, Dh, A, refpropm):
    # RPI = RefPropInterface()
    # refpropm = RPI.refpropm
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

    # From the book reference
    # HEAT TRANSFER AND AND FLUID FLOW IN MINICHANNELS AND MICROCHANNELS
    # Kandlikar et al, Elsevier p. 283
    # Martinelli parameter
    X = ((1-VQ)/VQ)**(0.9)*(Dvap/Dliq)**(0.5)*(Vvap/Vliq)**0.1
    # Modified Chisholm correlation for high flows (6.49)
    # Data from Wang et al (1997)
    # Vapor pressure multiplier
    TPVap2 = 1+9.4*(X**(0.62))+0.564*(X**(2.45))
    if MFLX < 200:
        C = 4.566*(X**(0.128))*(REliqO**(0.938))*((Dvap/Dliq)**(2.15))*((Vliq/Vvap)**(5.1))
        TPVap2 = 1+C*X+C*X*X
    dP=TPVap2*dPvap
    
#    dP=TPmult2*dPliqO #2-phase pressure gradient bar/m

    e=F_Void_8(MFLX,VQ,Dliq,Dvap,Isft) #Void fraction
    rm=(Dvap*e+Dliq*(1-e))*A #Relative mass kg/m

    return dP, rm
