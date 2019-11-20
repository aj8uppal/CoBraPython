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
  Thermal diffusivity [cm^2/s]
^   Prandtl number [-]'''

class RefPropInterface:
    def __init__(self):
        RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
        RP.SETPATHdll(os.environ['RPPREFIX'])
        self.UNIT = 0#RP.GETENUMdll(0,"MKS").iEnum
        self.RP = RP
    def refpropm(self, output, iName1, iValue1, iName2, iValue2, fluid, unit=None):
        return self.RP.REFPROPdll(fluid, iName1+iName2, output, unit or self.UNIT, 0, 0, iValue1, iValue2, [1.0]).Output[0]
