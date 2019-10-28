#VERIFIED

def MolarMass(fluid):
    #Commented until refprop is fixed, arbitrary value for now
    # Dco2=refpropm('D','T',293.15,'P',1e2,'CO2');
    # D=refpropm('D','T',293.15,'P',1e2,fluid);
    Dco2 = 1.98
    D = 0.807

    mm = 44.01*D/Dco2
    #R134a: 102.03 kg/kmol
    #CO2:44.01 kg/kmol
    #R218: 188.02 kg/kmol
    return mm
