import pandas as pd
import numpy as np
from helper import roundup

class CoolingBranch_v1a:
    def __init__(self, Fluid, Tsp, Tsc, Tsh, MF, HF, config, dL, name, plt, HFf):
        #config = {filename: , sheetname: }
        self.Fluid = Fluid
        self.Tsp = Tsp
        self.Tsc = Tsc
        self.Tsh = Tsh
        self.MF = MF
        self.HF = HF
        self.dL = dL
        self.name = name
        self.plt = plt
        self.HFf = HFf
        df = pd.read_excel(config['filename'], sheet_name=config['sheetname'])
        getCol = lambda N: df[df.columns[N if N < len(df.columns) else 0]]
        N_D=2;
        N_A=3;
        N_L=4;
        N_An=5;
        N_Ep=6;
        N_Q=7;
        N_HXnode=8;
        N_HXflowdir=9;
        N_HXcond=10  ;
        N_ISOcond=11;
        N_Tenv=12;
        N_HF=13;
        N_QPower=18;
        self.D=getCol(N_D)/1e3 #extract diameters
        self.A=getCol(N_A)/1e6 #extract flow crossectional areas
        self.L=getCol(N_L) #extract lengths
        self.An=getCol(N_An) #extract inclinations
        self.Ep=getCol(N_Ep)/1e6 #extract  roughnesses
        self.Q=getCol(N_Q) #extract heater labels
        self.HXnode=getCol(N_HXnode) #extract heat exchange nodes
        self.HXflowdir=getCol(N_HXflowdir) #extract heat exchanger flow direction
        self.HXcond=getCol(N_HXcond) #extract heat exchanger conductance
        self.ISOcond=getCol(N_ISOcond) #extract insulation conductance
        self.Tenv=getCol(N_Tenv) #extract Tenv
        self.HFinp=getCol(N_HF)
        self.QPower=getCol(N_QPower)
    def initialize_arrays(self):
        self.HFLX = [None]*(len(self.Q))
        for x in range(0, len(self.Q)):
            if self.Q[x] == 0:
                self.HFLX[x] = 0
            elif np.isnan(self.HF):
                self.HFLX[x]=self.HFinp[x]/(self.L[x]*np.pi*self.D[x])
            elif self.Q[x] == -1:
                self.HFLX[x]=self.QPower[x]/(self.L[x]*np.pi*self.D[x]) #Heatflux array
            else:
                self.HFLX[x]=self.HF[self.Q[x]-1]/(self.L[x]*np.pi*self.D[x]) #Heatflux array
    def fine_config(self):
        self.Lf=[0]
        self.SP=[1]
        self.Df=[0]
        self.Anf=[0]
        self.HFLXf=[0]
        self.HXnode_f=[0]
        self.HXflowdir_f=[0]
        self.HXcond_f=[0]
        self.HFLXf_app = []
        self.ISOcond_f=[0]
        self.Tenv_f=[0]
        self.Af_per = []
        self.Epf = []
        for x in range(len(self.L)):
            if np.isnan(self.HXnode[x]):
                N = self.L[x]/self.dL
            else:
                N = self.L[self.HXnode[x]]/self.dL
            N=roundup(N)
            dLtm=self.L[x]/N
            for y in range(len(self.Lf), int(N)+len(self.Lf)):
                self.Lf.append(self.Lf[-1]+dLtm)
                self.Df.append(self.D[x])
                self.Af_per.append(np.pi*self.Df[-1]*dLtm)
                self.Anf.append(self.An[x])
                self.Epf.append(self.Ep[x])
                self.HFLXf_app.append(self.HFLX[x])
                self.HXnode_f.append(self.HXnode[x])
                self.HXflowdir_f.append(self.HXflowdir[x])
                self.HXcond_f.append(self.HXcond[x])
                self.ISOcond_f.append(self.ISOcond[x])
                self.Tenv_f.append(self.Tenv[x])
            self.SP.append(len(self.Lf))
    def redefine(self):
        self.HXnode_fCopy=self.HXnode_f #Redefine fine HX node arrangement taking flowdirection into account.
        for x in range(1, len(self.SP)):
            z=1;
            if self.HXflowdir_f[self.SP[x]] == -1:  #counter current
                N=self.HXnode_f[self.SP[x]-1];
                for y in range(self.SP[x], self.SP[x-1]+1, -1):
                   self.HXnode_f[y-1]=self.SP[N-1]+z;
                   z=z+1;
            elif self.HXflowdir_f[self.SP[x]] == 1: #co current
                N=self.HXnode_f[self.SP[x]];
                for y in range(self.SP[x-1]+1, self.SP[x]):
                   self.HXnode_f[y-1]=self.SP[N]+z;
                   z=z+1;
        leng=1;
        if not leng == len(self.Lf):
            self.dP = []
            self.HTC = []
            for x in range(len(self.LF)): # set initial arrays for itteration start
                self.dP.append(1);
                self.HTC.append(1);
            # P(x)=refpropm('P','T',Tsp+273.15,'Q',0,Fluid)*1e-2;
        #     if strcmp(Tsc,'sat')==1 || strcmp(Tsc(1),'x')
        #     T(x)=Tsp;
        #     else
        #     T(x)=Tsp+Tsc;
        #     end
        #     H(x)=refpropm('H','T',Tsp+273.15,'Q',0,Fluid);
        #     MFLXf(x)=MF/(0.25*pi*self.Df(x)^2);
        #
        # end
