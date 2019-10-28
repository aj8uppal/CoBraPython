def ThomeCorrelation_Con(Fluid, P, H, MFLX, HFLX, dh, A, Ph):
  fluid = fluid
    G=MFLX #kg/m2s
    d=Dh #Equivalent diameter.
    p=P*1e2
    q=abs(HFLX) #w/m2K

    #Commented until refprop is fixed, arbitrary value for now
    #Inlet Vapor Quality
    # x=refpropm('Q','P',p,'H',H,fluid);
    #
    # #Enthalpy [J/kg]
    #
    # Hl=refpropm('H','P',p,'Q',0,fluid);
    # Hv=refpropm('H','P',p,'Q',1,fluid);
    #
    # #Density [kg/m3]
    # Dl=refpropm('D','P',p,'Q',0,fluid);
    # Dv=refpropm('D','P',p,'Q',1,fluid);
    #
    # #Viscosity [Pa.s]
    # Vl=refpropm('V','P',p,'Q',0,fluid);
    # Vv=refpropm('V','P',p,'Q',1,fluid);
    #
    # #Thermal Conductivity [W/m.K]
    # Kl=refpropm('L','P',p,'Q',0,fluid);
    # Kv=refpropm('L','P',p,'Q',1,fluid);
    #
    # #Specific Heat [J/kg.K]
    # CPl=refpropm('C','P',p,'Q',0,fluid);
    # CPv=refpropm('C','P',p,'Q',1,fluid);
    #
    # #Surface Tension [N/m]
    # ST=refpropm('I','P',p,'Q',0,fluid);

    x = Hl = Hv = Dl = Dv = Vl = Vv = Kl = Kv = CPl = CPv = ST = 2

    #Intermittent to Annular Flow Transition Boundary
    xia= F_xia_18con(Dl,Dv,Vl,Vv);
    eia=F_Void_8(G,xia,Dl,Dv,ST);

    #Void fraction
    e=F_Void_8(G,x,Dl,Dv,ST);

    #Get Transition Boundary - Stratified to Stratified-Wavy (S-SW)
    Gstrat=F_Gstrat_17con(d,A,x,e,Dl,Dv,Vl);


    #Get Transition Boundary - Stratified-Wavy to Intermittent and Annular (SW-I/A)
    Gwavy=F_Gwavy_16con(G,d,A,x,e,Dl,Dv,ST);

    #Transition Boundary - Intermittent to Bubbly Flow (I-B)
    Gbub=F_Gbub_22con(d,A,x,e,Dl,Dv,Vl);

    #Transition Boundary - Annular and intermittent flow to Mist Flow (AI-M)
    Gmist=F_Gmist_19con(G,d,A,x,e,Dl,Dv,ST);

function[dP,HTC,x,rm,flowpattern] = ThomeCorrelation_Con(Fluid,P,H,MFLX,HFLX,Dh,A,Ph)





%% flow definition

%Fully Stratified Flow
if G<=Gstrat
flowpattern='Strat';
FilmAngle=F_StratAngle_13(e);

%Stratified-Wavy Flow (SW)
elseif G>=Gstrat && G<=Gwavy
flowpattern='SW';
FilmAngle=F_FilmAngle_8130con( G,d,A,x,e,Dl,Dv,Vl,ST);

%Intermittent Flow (I)
elseif G>=Gwavy && G<=Gmist && x<=xia && G<=Gbub
flowpattern='Int';
FilmAngle=0;

%Bubbly flow
elseif  G<=Gmist && x<=xia && G>Gbub
flowpattern='Bub';
FilmAngle=0;

%Annular Flow (A)
elseif G>=Gwavy && G<=Gmist && x>xia
flowpattern='Annu';
FilmAngle=0;

%Mist Flow (M)
elseif  G>Gmist && G>=Gstrat
flowpattern='Mist';
FilmAngle=0;


else flowpattern='NotIdentified';
end
%%
%%
% if x<=1
if strcmp(flowpattern,'Strat')==1
    DPf=F_DPstrat_52(G,d,x,e,Dl,Dv,Vl,Vv,ST );
    htp=F_HTC_8123con( G,q,d,A,x,e,FilmAngle,Dl,Dv,CPl,Kl,Hl,Hv,Vl,ST);

elseif strcmp(flowpattern,'SW')==1
    DPf=F_DPsw_39(G,d,A,x,e,Dl,Dv,Vl,Vv,ST);
    htp=F_HTC_8123con( G,q,d,A,x,e,FilmAngle,Dl,Dv,CPl,Kl,Hl,Hv,Vl,ST);

elseif strcmp(flowpattern,'Int')==1
    DPf=F_DPslug_int_35(G,d,x,e,Dl,Dv,Vl,Vv,ST);
     htp=F_HTC_8123con( G,q,d,A,x,e,FilmAngle,Dl,Dv,CPl,Kl,Hl,Hv,Vl,ST);

elseif strcmp(flowpattern,'Annu')==1
    DPf=F_DPannu_29(G,d,x,e,Dl,Dv,Vv,ST);
    htp=F_HTC_8123con( G,q,d,A,x,e,FilmAngle,Dl,Dv,CPl,Kl,Hl,Hv,Vl,ST);

elseif strcmp(flowpattern,'Mist')==1
    DPf=F_DPmist_45(G,d,x,Dl,Dv,Vl,Vv);
    htp=F_HTC_8123con( G,q,d,A,x,e,FilmAngle,Dl,Dv,CPl,Kl,Hl,Hv,Vl,ST);

elseif strcmp(flowpattern,'Bub')==1
    DPf=F_DPbub_56(G,d,x,e,Dl,Dv,Vl,Vv,ST );
    htp=F_HTC_8123con( G,q,d,A,x,e,FilmAngle,Dl,Dv,CPl,Kl,Hl,Hv,Vl,ST);

elseif strcmp(flowpattern,'NotIdentified')==1
    DPf=0;
    htp=0;
end
DPmom=F_DPmomentum_28(G,q,A,Ph,x,e,Dl,Dv,Hl,Hv,ST);
dP=DPf+DPmom;
%%

dP=1e-5*dP; %bar
HTC=real(htp);

rm=(Dv*e+Dl*(1-e))*A; %Relative mass kg/m



end
