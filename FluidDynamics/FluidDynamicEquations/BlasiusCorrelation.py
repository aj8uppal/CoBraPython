#VERIFIED

def BlasiusCorrelation(Re):
    '''Calculation of the Blasius friction Factors for the Darcy Weisbach formula for pressure drop. Input Reynolds Re'''
    if Re < 2300:
        Fd = 64/Re
        State = 'Lam'
    elif Re < 20000:
        Fd=0.3164/Re**0.25
        State = 'Tur'
    else:
        Fd=0.184/Re**0.2
        State = 'Tur'
    return Fd, State
