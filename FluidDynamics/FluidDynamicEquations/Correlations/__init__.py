from cmath import log, sqrt


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
