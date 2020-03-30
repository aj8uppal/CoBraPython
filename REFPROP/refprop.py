import os, numpy as np
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary

'''0   Refprop DLL version number
A   Speed of sound [m/s]
B   Volumetric expansivity (beta) [1/K]
C   Cp [J/(kg K)]
D   Density [kg/m3]
E   dP/dT (along the saturation line) [kPa/K]
F   Fugacity [kPa] (returned as an array)
G   Gross heating value [J/kg]
H   Enthalpy [J/kg]
I   Surface tension [N/m]
J   Isenthalpic Joule-Thompson coeff [K/kPa]
K   Ratio of specific heats (Cp/Cv) [-]
L   Thermal conductivity [W/(m K)]
M   Molar mass [g/mol]
N   Net heating value [J/kg]
O   Cv [J/(kg K)]
P   Pressure [kPa]
Q   Quality (vapor fraction) (kg/kg)
R   d(rho)/dP (constant T) [kg/kPa]
S   Entropy [J/(kg K)]
T   Temperature [K]
U   Internal energy [J/kg]
V   Dynamic viscosity [Pa*s]
W   d(rho)/dT (constant p)[kg/(m^3 K)]
X   Liquid phase & gas phase comp.(mass frac.)
Z   Compressibility factor
+   Liquid density of equilibrium phase
-   Vapor density of equilibrium phase
!   dH/d(rho) (constant T) [(J/kg)/(kg/m^3)]
@   dH/dT (constant rho) [J/(kg K)]
#   dP/dT (constant rho) [kPa/K]
$   Kinematic viscosity [cm^2/s]
%   Thermal diffusivity [cm^2/s]
^   Prandtl number [-]'''

class RefPropInterface:
    def __init__(self):
        RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
        RP.SETPATHdll(os.environ['RPPREFIX'])
        self.UNIT = 7 #RP.GETENUMdll(0,"MKS").iEnum
        self.RP = RP
    def changeUnits(self, inputUnit, factor=False):
        outputUnit = inputUnit
        #CritP -->
        if inputUnit == 'V':
            outputUnit = 'VIS'
            factor = factor if not factor else 1e-6
        elif inputUnit == 'C':
            outputUnit = 'CP'
            factor = factor if not factor else 1e3
        elif inputUnit == 'I':
            outputUnit = 'STN'
            factor = factor if not factor else 1e-3
        elif inputUnit == 'L':
            outputUnit = 'TCX'
        elif inputUnit == 'H':
            factor = factor if not factor else 1e3
        return outputUnit, 1 if not factor else factor
    def refpropm(self, output, iName1, iValue1, iName2, iValue2, fluid, unit=None):
        factor = 1
        output, factor = self.changeUnits(output,factor=True)
        iName1, _ = self.changeUnits(iName1)
        iName2, _ = self.changeUnits(iName2)
        iValue1 = iValue1*1e-3 if iName1 == 'H' else iValue1
        iValue2 = iValue2*1e-3 if iName2 == 'H' else iValue2
        if output == 'Q':
            retval = self.RP.REFPROPdll(fluid, iName1+iName2, 'QMASS', 7, 0, 0, iValue1, iValue2, [1.0]).q
            if retval < 1e-4:
                return 1e-4
            elif retval > (1-1e-4):
                return (1-1e-4)
            else:
                return retval
        return factor*self.RP.REFPROPdll(fluid, iName1+iName2, output, unit or self.UNIT, 0, 0, iValue1, iValue2, [1.0]).Output[0]
    def closeTo(self, a1, a2, epsilon=0.01):
        return (abs(a1-a2) < abs((a1+a2)/2)*epsilon)
    def test(self):
        #complete: I, C, L, V, D, H, P, Q
        assert(self.closeTo(self.refpropm('I','P',1000,'Q',0,'CO2'), 0.0127))
        assert(self.closeTo(self.refpropm('C','P',1000,'Q',0,'CO2'), 2.0111e+03))
        assert(self.closeTo(self.refpropm('L','P',1000,'Q',0,'CO2'), 0.1595, epsilon=0.02))
        assert(self.closeTo(self.refpropm('V','P',1000,'Q',0,'CO2'), 1.9415e-04))
        assert(self.closeTo(self.refpropm('D','P',1000,'Q',0,'CO2'), 1.1169e+03))
        assert(self.closeTo(self.refpropm('H','P',1000,'Q',0,'CO2'), 1.1266e+05))
        assert(self.closeTo(self.refpropm('Q','P',1000,'H',2e5,'CO2'), 0.2707))
        assert(self.closeTo(self.refpropm('D','P',1000,'H',2e5,'CO2'), 90.3934))
        assert(self.closeTo(self.refpropm('V','P',1000,'H',2e6,'CO2'), 5.5715e-05))
        assert(self.closeTo(self.refpropm('L','P',1000,'H',2e6,'CO2'), 0.1040, epsilon=0.025))
