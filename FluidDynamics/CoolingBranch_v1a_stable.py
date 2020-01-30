import pandas as pd
import numpy as np
from helper import roundup
from math import pi
from inspect import signature
import sys
sys.path.append('../REFPROP')
from refprop import RefPropInterface
sys.path.append('PressureDropAndHeatTransfer')
from time import time
from dPandHTC import *
import os

#use np arrays

import warnings
warnings.filterwarnings('ignore')

class CoolingBranch_v1a:
    refpropm = RefPropInterface().refpropm
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
        self.nargin = len(signature(CoolingBranch_v1a).parameters)
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
        # print(self.ISOcond)
        self.Tenv=getCol(N_Tenv) #extract Tenv
        self.HFinp=getCol(N_HF)
        self.QPower=getCol(N_QPower)
    def start(self):
        self.initialize_arrays()
        self.fine_config()
        self.redefine()
    def initialize_arrays(self):
        self.HFLX = [float('nan')]*(len(self.Q))
        for x in range(0, len(self.Q)):
            if self.Q[x] == 0:
                self.HFLX[x] = 0
            elif all([np.isnan(i) for i in self.HF]):
                self.HFLX[x]=self.HFinp[x]/(self.L[x]*np.pi*self.D[x]) #r_tube
            elif self.Q[x] == -1:
                self.HFLX[x]=self.QPower[x]/(self.L[x]*np.pi*self.D[x]) #Heatflux array
            else:
                self.HFLX[x]=self.HF[self.Q[x]-1]/(self.L[x]*np.pi*self.D[x]) #Heatflux array
    def fine_config(self):
        self.Lf=[0]
        self.SP=[1]
        self.Df=[float('nan')]
        self.Anf=[float('nan')]
        self.HFLXf=[float('nan')]
        self.HXnode_f=[float('nan')]
        self.HXflowdir_f=[float('nan')]
        self.HXcond_f=[float('nan')]
        self.HFLXf_app = [0]
        self.ISOcond_f=[float('nan')]
        self.Tenv_f=[float('nan')]
        self.Af_per = [0]
        self.Epf = [0]
        for x in range(len(self.L)):
            if np.isnan(self.HXnode[x]):
                N = self.L[x]/self.dL
            else:
                N = self.L[self.HXnode[x]]/self.dL
            N=roundup(N)
            dLtm=self.L[x]/N
            for y in range(len(self.Lf), int(N)+len(self.Lf)):
                self.Lf.append(self.Lf[y-1]+dLtm)
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
        self.VQ = [float('nan')]*len(self.Lf)
        self.rm = [float('nan')]*len(self.Lf)
        self.dH = [float('nan')]*len(self.Lf)
        self.Tw = [float('nan')]*len(self.Lf)
        self.HFLXf = [float('nan')]*len(self.Lf)
        self.State = ['']*len(self.Lf)
    def redefine(self):
        self.HXnode_fCopy=self.HXnode_f #Redefine fine HX node arrangement taking flowdirection into account.
        for x in range(1, len(self.SP)):
            z=1;
            # print(self.SP, x, self.HXflowdir_f)
            if self.HXflowdir_f[self.SP[x]-1] == -1:  #counter current
                N=self.HXnode_f[self.SP[x]-1];
                for y in range(self.SP[x], self.SP[x-1]+1, -1):
                   self.HXnode_f[y-1]=self.SP[N-1]+z;
                   z=z+1;
            elif self.HXflowdir_f[self.SP[x-1]] == 1: #co current
                N=self.HXnode_f[self.SP[x]];
                for y in range(self.SP[x-1]+1, self.SP[x]):
                   self.HXnode_f[y-1]=self.SP[N-1]+z;
                   z=z+1;
        leng=1;
        if not leng == len(self.Lf):
            self.dP = []
            self.HTC = []
            self.P = []
            self.T = []
            self.H = []
            self.MFLXf = []
            for x in range(len(self.Lf)): # set initial arrays for itteration start
                self.dP.append(1);
                self.HTC.append(1);
                self.P.append(self.refpropm('P','T',self.Tsp+273.15,'Q',0,self.Fluid)*1e-2);
                if self.Tsc == 'sat':
                    self.T.append(self.Tsp)
                else:
                    self.T.append(self.Tsp+self.Tsc)
                self.H.append(self.refpropm('H','T',self.Tsp+273.15,'Q',0,self.Fluid))
                self.MFLXf.append(self.Df[x] if not self.Df[x] else self.MF/(0.25*pi*self.Df[x]**2));

    def run(self):
        itt=0;
        converge=10000;
        HXstart=0;
        convlimit=40;
        conv_repeat=0;
        conv_repeat_limit=3;
        self.HFLXf_hx = [float('nan')]*len(self.Lf)
        self.HFf_hx = [0]*len(self.Lf)
        self.HFLXf_env = [0]*len(self.Lf)
        self.Hconv = []
        self.HTCconv = []
        ittstop = 400
        while (abs(converge)>convlimit or conv_repeat<conv_repeat_limit+1) and itt<ittstop:
            foo = True
            itt+=1
            print("Iteration Round: {}".format(itt))
            print("Iteration offset: {} (Stops at {})".format(converge, convlimit))
            print("Enthalpy: {} J/kg".format(max(self.H)))
            print("Iteration offset: {} (Stops at {})".format(conv_repeat, conv_repeat_limit))
            Pprev = self.P
            Tprev = self.T
            HTCprev = self.HTC
            Hprev = self.H
            if converge < 1000:
                HXstart = 1
            times = {1: 0, 2: 0, 3: 0, 4: 0}
            startTime = time()
            # rows, columns = os.popen('stty size', 'r').read().split()
            columns, rows = os.get_terminal_size()
            for x in range(1, len(self.Lf)): #2:length(Lf)
                elapsed = time()-startTime
                total = elapsed*len(self.Lf)/x
                remaining = total-elapsed
                prefix = "\r{}/{} ({}%) [Elapsed time: {} Remaining time: {}] ".format(x+1, len(self.Lf), round(x/len(self.Lf)*100), round(elapsed), round(remaining))
                progressbar_len = int(columns) - (len(prefix) + 2)
                progressbar = "".join(["#" if j <= x*progressbar_len//len(self.Lf) else "." for j in range(progressbar_len)])
                sys.stdout.write(prefix+"["+progressbar+"]");
                sys.stdout.flush()
                # print(len(self.Lf))
                # print("{}/{}".format(x, len(self.Lf)))
                #section 1
                # start = time()
                if np.isnan(self.ISOcond_f[x]): #find environmental leak
                    self.HFLXf_env[x] = 0
                else:
                    HFLXf_env_prev = self.HFLXf_env[x]
                    self.HFLXf_env[x] = ((self.Tenv_f[x]-self.T[x])*(((self.HTC[x]*self.Af_per[x])**(-1)+(self.ISOcond_f[x]*(self.Lf[x]-self.Lf[x-1]))**(-1))**(-1)))/self.Af_per[x]
                    n = 0.1 #average heat exchange load over 1/n samples (To stabilze conversion)
                    self.HFLXf_env[x] = HFLXf_env_prev*(1-n)+self.HFLXf_env[x]*n
                # times[1]+=(time()-start)
                #section 2
                # start = time()
                if np.isnan(self.HXcond_f[x]):
                    if np.isnan(self.HFLXf_hx[x]) or abs(self.HFLXf_hx[x]) < 0:
                        self.HFLXf_hx[x] = 0
                elif HXstart == 1: #if heat exchange node exist calculate heat exchange by temperature difference
                    HFf_hx_prev = HFf_hx[x]
                    self.HFf_hx[x]=((self.T[self.HXnode_f[x]-1]-self.T[x])*(((self.HTC[x]*self.Af_per[x])**(-1)+(self.HXcond_f[x]*(self.Lf[x]-self.Lf[x-1]))**(-1)+(self.HTC[self.HXnode_f[x]-1]*self.Af_per[self.HXnode_f[x]-1])**(-1))**(-1)))
                    n = 0.1 #average heat exchange load
                    self.HFf_hx[x] = HFf_hx_prev*(1-n)+self.HFf_hx[x]*n
                    self.HFf_hx[self.HXnode_f[x]-1] = -HFf_hx[x]
                    self.HFLXf_hx[x]=self.HFf_hx[x]/self.Af_per[x]
                    self.HFLXf_hx[HXnode_f[x]-1]=self.HFf_hx[HXnode_f[x]-1]/self.Af_per[self.HXnode_f[x]-1]
                else:
                    HFLXf_hx[x] = 0
                # times[2]+=(time()-start)
                #section 3
                # start = time()
                if self.nargin >= 11 and len(self.HFf) == len(self.Lf):
                    self.HFLXf[x] = self.HFf[x]/self.Af_per[x]+self.HFLXf_hx[x]
                else:
                    self.HFLXf[x]=self.HFLXf_env[x]+self.HFLXf_hx[x]+self.HFLXf_app[x]
                    if all(np.isnan(self.HF)):
                        self.HFf = [sum(self.HFinp)]
                    else:
                        self.HFf = [sum(self.HF)]
                # times[3]+=(time()-start)
                #section 4
                # start = time()
                if self.Fluid == 'CO2' and self.H[x] < 9e4: #avoid enthalpy to get into freezing area for cO2.
                    self.H[x] = 9e4
                # print("Entering dPandHTC...")
                newDP, newHTC, newVQ, newRM, newState, newT = dPandHTC(self.Fluid, self.P[x], self.H[x], self.MFLXf[x], self.HFLXf[x], self.Df[x], 0.25*pi*self.Df[x]**2, pi*self.Df[x], self.Epf[x], self.Anf[x], self.Tsh, self.refpropm);
                # print("{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}".format(self.Fluid, self.P[x], self.H[x], self.MFLXf[x], self.HFLXf[x], self.Df[x], 0.25*pi*self.Df[x]**2, pi*self.Df[x], self.Epf[x], self.Anf[x], self.Tsh))
                # times[4]+=(time()-start)
                # print("Exiting dPandHTC...")
                self.dP[x] = newDP
                self.HTC[x] = newHTC
                self.VQ[x] = newVQ
                self.rm[x] = newRM
                self.State[x] = newState
                self.T[x] = newT
                self.dH[x] = (self.HFLXf[x]*pi*self.Df[x]*(self.Lf[x]-self.Lf[x-1]))/self.MF #calculate enthalpy difference, mass not mol
                # print(self.dH[x
                if(abs(self.dH[x]) > 1e6):
                    foo = False
                    break
                self.Tw[x] = self.T[x]+(self.HFLXf[x]/self.HTC[x]) #wall temperature
                self.H[x] = self.H[x-1]+self.dH[x] #calculate new enthalpy
                # print(self.dH[x])

            # print(times)
                #end for loop
            # foo = max([i for i in self.dH if not np.isnan(i)])
            # if foo > 1e6:
            #     break
            if not foo:
                break
            self.P[-1] = self.refpropm('P','T',self.Tsp+273.15,'Q',1,self.Fluid)*1e-2;
            for x in range(len(self.Lf)-2, -1, -1):
                self.P[x] = self.P[x+1]+self.dP[x+1]*(self.Lf[x+1]-self.Lf[x]) #calculate new pressure
            if self.Tsc == 'sat': #force inlet to be saturated
                self.H[0] = self.refpropm('H','P',min(73, self.P[0])*1e2,'Q',0,self.Fluid) #saturated inlet enthalpy
                # if self.P[0] > 73:
                #     self.H[0] = refpropm('H','P',73*1e2,'Q',0,self.Fluid)
                # else:
                #     self.H[0] = refpropm('H','P',self.P[0]*1e2,'Q',0,self.Fluid) #saturated inlet enthalpy
                self.VQ[0] = 0
                self.T[0] = float('nan')
            elif type(self.Tsc) == str and self.Tsc[0] == 'x':
                xi = int(self.Tsc[1:])/100
                self.H[0] = self.refpropm('H','P',self.P[0]*1e2,'Q',xi,self.Fluid) #calculate saturated inlet enthalpy
                self.VQ[0] = xi
                self.T[0] = float('nan')
            else:
                self.H[0] = self.refpropm('H','T',self.T[0]+273.15,'P',self.P[0]*1e2,self.Fluid) #calculate inlet enthalpy based on given inlet temperature

            self.Tw[0] = float('nan')
            # Pconv=P-Pprev;
            # Tconv=T-Tprev;
            self.Hconv = np.array(self.H)-np.array(Hprev)
            self.HTCconv = np.array(self.HTC)-np.array(HTCprev)
            converge_prev = converge
            converge = max(self.Hconv) #use enthalpy to converge
            if abs(converge) < convlimit:
                conv_repeat+=1
            else:
                conv_repeat = 0
            #end while loop
    def plot(self):
        pass

x = CoolingBranch_v1a('CO2', -25, 0, 0, 5*1.516e-3, [134.4,100], {'filename': '../CoBraV1a_example.xlsx', 'sheetname': 'Example'}, 10e-3, "Plot", 0, [3600])
x.start()
x.run()
print("Completed.")

