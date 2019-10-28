from cmath import log, sqrt
from numpy import pi


def F_Void_8(G, x, Dl, Dv, ST):
    g = 9.81
    # Void Fraction (Rouhani-Axelsson drift flux model)
    #e=(x/Dv)*((1+0.12*(1-x))*((x/Dv)+((1-x)/Dl))+((1.18*(1-x)*(g*ST*(Dl-Dv))^(1/4))/(G*Dl^(1/2))))^(-1)

    # Homogeneous Void Fraction
    eH=1/(1+((1-x)/x)*(Dv/Dl));

    #Rouhani Void Fraction
    er=(x/Dv)*((1+0.12*(1-x))*((x/Dv)+((1-x)/Dl))+((1.18*(1-x)*(g*ST*(Dl-Dv))**(1/4))/(G*Dl**(1/2))))**(-1)

    #Logarithmic Mean Void Fraction
    e = (eH-er)/log(eH/er)
    e = er
    return e

def F_xia_18con(Dl, Dv, Vl, Vv):
    xia=(0.2914*(Dv/Dl)**(-1/1.75)*(Vl/Vv)**(-1/7)+1)**(-1) #0.2914 replaced 1.8^(1/0.875)
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
    # Gstrat=((226.3^2*Ald_xia*Avd_xia^2*Dv*(Dl-Dv)*Vl*g)/(xia^2*(1-xia)*pi^3))^(1/3);
    # else
    Gstrat=((226.3**2*Ald*Avd**2*Dv*(Dl-Dv)*Vl*g)/(x**2*(1-x)*pi**3))**(1/3)+20*x;   #20*x is new
    return Gstrat
