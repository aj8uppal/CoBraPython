from numpy import log10

def ColebrookEquation(Dh, Ep, Re):
    Fd = 1e-3
    aa = 1
    while aa > 1e-5:
        bb = -2*log10(((Ep/Dh)/3.7)+2.51/(Re*Fd**0.5))
        Fdprev = Fd;
        Fd = bb**(-2);
        aa = abs(Fd-Fdprev);
    return Fd
