import sys
sys.path.append('../../REFPROP')
from refprop import RefPropInterface

def SteinerTaborek(fluid,P,H,MFLX,HFLX,Dh,A,Ph,g):
    '''Data book
    Chapter10 boilinig heat transfer inside plain tubes 10.3.4 steiner-taborek
    asymptotic model
    '''
    q=HFLX;
    z=5 #m       length of the tube?    7.64<Z<8.84m
    RPI = RefPropInterface()
    refpropm = RPI.refpropm
    G=MFLX;
    # Tsat=refpropm('T','P',P,'Q',0,fluid);
    Tsat = 2
    Tk=Tsat;
    Pcri=73.8;  # bar
    P1=P/100;   # bar
    Pr=P1/Pcri;

    #Enthalpy [J/kg]
    Hl=refpropm('H','P',P,'Q',0,fluid);
    Hv=refpropm('H','P',P,'Q',1,fluid);

    #Density [kg/m3]
    Dl=refpropm('D','P',P,'Q',0,fluid);
    Dv=refpropm('D','P',P,'Q',1,fluid);

    #Viscosity [Pa.s]
    Vl=refpropm('V','P',P,'Q',0,fluid);
    Vv=refpropm('V','P',P,'Q',1,fluid);


    #Thermal Conductivity [W/m.K]
    Kl=refpropm('L','P',P,'Q',0,fluid);
    Kv=refpropm('L','P',P,'Q',1,fluid);

    #Specific Heat [J/kg.K]
    CPl=refpropm('C','P',P,'Q',0,fluid);
    CPv=refpropm('C','P',P,'Q',1,fluid);

    #Surface Tension [N/m]
    ST=refpropm('I','P',P,'Q',0,fluid);
    Hl = Hv = Dl = Dv = Vl = Vv = Kl = Kv = CPl = CPv = ST = 2
    Hv = 2*Hl

    xcrit=0.50;
    #Chapter 18.6 Critical heat flux in vertical channel
    # xcrit=x_critver(fluid,P,H,MFLX,Dh,z);
    # fprintf('xcrit is %f\n',xcrit);

    ro=0.3*10e-6;
    HTC_xt=[F_HTC_xtver(G,Dh,Vl,CPl,Kl)]
    q_ONB=F_q_ONBver(Dv,Tk,ST,HTC_xt,Hl,Hv,ro)

    HTC_nbo=18890;#for CO2, with property at Tcrit
    qo=150000;
    M=44.01;
    Fnb=[F_Fnbver(Pr,Dh,q,qo,M)]

     # x=refpropm('Q','P',P,'H',H,fluid);
     x=2
    #  e=F_voidver(x,Dv,g,Dh,Dl,G,ST);           %  Data book Chapter 17 Page17-18 Rouhani-Axelsson Correlation

    if x==0 && q<q_ONB:
        HTC_xt=[F_HTC_xtver(G,Dh,Vl,CPl,Kl)]
        HTC_tp=HTC_xt;
    elif x==0 && q>q_ONB:
        Ftp=((1-x)**1.5+1.9*x**0.6*(Dl/Dv)**0.35)**1.1;
        [HTC_xt]=F_HTC_xtver(G,Dh,Vl,CPl,Kl);
        HTC_tp=((HTC_nbo*Fnb)**3+(HTC_xt*Ftp)**3)**(1/3);
    elif x<xcrit && q>q_ONB:
        Ftp=((1-x)**1.5+1.9*x**0.6*(Dl/Dv)**0.35)**1.1;
        [HTC_xt]=F_HTC_xtver(G,Dh,Vl,CPl,Kl);
        HTC_tp=((HTC_nbo*Fnb)**3+(HTC_xt*Ftp)**3)**(1/3);
    elif x>=xcrit && x<1:   # mist flow
        HTC_mist,x=F_mistver(fluid,P,x,G,Dh) #mist flow HTC by 18.5.2 Groeneveld Method Page18-9
        HTC_tp= HTC_mist;
    else:
        #x==1 && q<q_ONB:
        HTC_Lt=[F_HTC_xtver(G,Dh,Vl,CPl,Kl)]
        HTC_Gt=[F_HTC_xtver(G,Dh,Vv,CPv,Kv)]
        Ftp=(((1-x)**1.5+1.9*x**0.6*(1-x)**0.01*(Dl/Dv)**0.35)**(-2.2)+(HTC_Lt/HTC_Gt*x**0.01*(1+8*(1-x)**0.7)*(Dl/Dv)**0.67)**(-2))**(-0.5);
        HTC_xt=[F_HTC_xtver(G,Dh,Vv,CPv,Kv)]
        HTC_tp=HTC_xt*Ftp;
    VQ=x;
    return HTC_tp,VQ
