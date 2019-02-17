%%%%%%%%             To understand the code read article        %%%%%%%%%%%%
%%  PARAMETER IDENTIFICATION PROBLEM IN THE HODGKIN AND HUXLEY MODEL      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 ESTIMATE CONDUCTANCES G_Na, G_K, G_L                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Model of Hodgkin and Huxley Temporal                    %
%             Euler Explicit:  To estimate  G_Na, G_K, G_L                 %
%                                                                          %
%|--------------------------------                                         %
%| CV_t=Iext-G_Na.m^a.h^b.(V-E_Na)-G_K.n^c.(V-E_K)-G_L.(V-E_L)             %
%| m_t=(1-m).AlphaNm(V)-m.BetaNm(V)                                        %
%| h_t=(1-h).AlphaNh(V)-h.BetaNh(V)                                        %
%| n_t=(1-n).AlphaNn(V)-n.BetaNn(V)                                        %
%| V(0)=V_0;   m(0)=m_0    h(0)=h_0   n(0)=n_0                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clc;close all;clear all                                                     
global a b c Ena Ek El nt dt m n h V Vp

%%-- Start: defining Ordinari Diferential Equation (ODE (7) of paper)   --%%
                                                                           %|
%----    Set the time interval  ([0 T] in the paper) ti=0 and tf=T      ---%|
ti=0;     tf=10;   nt=500;                                                  %|
t=linspace(ti,tf,nt) ;  dt = t(2)-t(1);                                    %|
                                                                           %|
%----                  known parameters                                 ---%|
C_M=1;  Iext=0.0;   Ena=115;   Ek=-12;   El=10.598;  a=3;   b=1;    c=4;        %|
                                                                           %|
%----                  initial conditions                               ---%|
V=zeros(1,nt);  m=zeros(1,nt); n=zeros(1,nt); h=zeros(1,nt);               %|
V(1)=-25;       m(1)=0.5;      n(1)=0.4;      h(1)=0.4;                    %|
                                                                           %|
%---                unknown parameters (Conductances)                      %|
G_Na=120;      G_K =36;        G_L=0.3;                                    %|
                                                                           %|
%---                guess initial (Conductances)                           %|
G_Nak=0;        G_Kk=0;        G_Lk=0;                                     %|
                                                                           %|
%---   perturbation of the electrical potential (in percentage )        ---%|
MaxErro=5/100;                                                             %|
%delta=10;                                                                 %|  
                                                                           %|
%---             for the stop criterion                                 ---%|
tau=1.01;                                                                   %|
%%---                              End                                ---%%%|

%Vp=importdata('F3-4.txt');                                                
%------------------------                    -----------------------------%|

%-----------            Calculating the exact Vexa              ------------%
Vexa=Vexata(V,G_Na,G_K,G_L,C_M,Iext);

%----------        Making the pertubation of Vexa in Vp         ------------%
Vp=Vexa + (-MaxErro+(2*MaxErro).*rand(1,nt)).*Vexa;

%---------  Calculing delta for the equation (9) of the paper   ------------%
delta=MaxErro*sqrt( dt^1*sum( (Vexa).^2 ) );



%%--------------------       k=========1       ------------%%
k=0;
                      
while(0==k || tau*delta<=ResiduoV(k)) 
%for  i=1:1
k=k+1;                                                                
%--------------    Calculing   Vk, Uk, mk nk and hk          ---------------%
  [Vk,Uk,mk,nk,hk]=Iaprox(G_Nak,G_Kk,G_Lk,C_M,Iext);

%-------------      Calculing of the residue: ||Vp-Vk||     ---------------%
  ResiduoV(k)=sqrt( dt^1*sum( (Vp-Vk).^2 ) );

%-------------        Calculing the Errors                  ---------------%
  Error(k)=sqrt( [(G_Nak-G_Na)^2+(G_Kk-G_K)^2+(G_Lk-G_L)^2 ]/...
               [(G_Na)^2+(G_K)^2+(G_L)^2] )*100;
  gna(k)  =G_Nak;  gk(k)   =G_Kk ;  gl(k)  =G_Lk;    
fprintf('%10.10f\t',k,G_Nak,G_Kk,G_Lk,Error(k),ResiduoV(k));  fprintf('\n\n\n');
%------               To the minimal error method           --------------%
  Wk=sum( dt*(Vp-Vk).^2 )/...
        ( ( sum( mk.^a.*hk.^b.*(Vk-Ena).*Uk ) )^2+...
        ( sum( nk.^c.*(Vk-Ek ).*Uk ) )^2+( sum( (Vk-El ).*Uk ) )^2 );

%-------------                   Print                      --------------%


%------------         Calculing the iteration k             --------------%
 G_Nak=G_Nak+Wk*dt*sum( mk.^a.*hk.^b.*(Vk-Ena).*Uk );%For G_Na
 G_Kk =G_Kk +Wk*dt*sum( nk.^c.*       (Vk-Ek ).*Uk );%For G_K
 G_Lk =G_Lk +Wk*dt*sum(               (Vk-El ).*Uk );%For G_L

end
fprintf('%10.10f\t',k,G_Nak,G_Kk,G_Lk,Error(k),ResiduoV(k));  fprintf('\n\n\n');
k=linspace(1,k,k);  ks=k(1:50:max(k)); if (max(k)~= max(ks) ) k= [ks max(k) ]; end;
gna=gna(k); gk=gk(k);  gl=gl(k); ResiduoV= ResiduoV(k); Error=Error(k);
save('Example-11.txt','Vexa','Vp','t','-ascii');
save('Example-12.txt','gna','gk','gl','ResiduoV','Error','k','-ascii');

%plot(t,Vp,t,Vexa,'r')
