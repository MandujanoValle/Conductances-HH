%---------------------------------------------------------------------------------------------------%
%---                                  Explicit Euler Method                                     ----%
%---------------------------------------------------------------------------------------------------%

function [V,U,m,n,h]=Iaprox(Gna,Gk,Gl,C,Iext)
global a b c Ena Ek El nt dt m  n  h V Vp

%-----                                Approximate solution of V                                -----%
d=dt/C;
for i=1:nt-1
AlphaMv(i)=( (25-V(i))/10 )/( exp( (25-V(i))/10 )-1 );  BetaMv(i) =4*exp(-V(i)/18);
AlphaNv(i)=(10-V(i))/( 100*( exp( (10-V(i))/10 )-1) );  BetaNv(i) =0.125*exp(-V(i)/80);
AlphaHv(i)=0.07*exp( -V(i)/20 );                        BetaHv(i) =1/( exp( (30-V(i))/10 )+1 );

V(i+1)  =V(i) + d* ( Iext - Gna*m(i)^a*h(i)^b*(V(i)-Ena) - Gk*n(i)^c*(V(i)-Ek) - Gl*(V(i)-El) );
m(i+1)  =m(i) + dt*( (1-m(i))*AlphaMv(i) - m(i)*BetaMv(i) );
n(i+1)  =n(i) + dt*( (1-n(i))*AlphaNv(i) - n(i)*BetaNv(i) );
h(i+1)  =h(i) + dt*( (1-h(i))*AlphaHv(i) - h(i)*BetaHv(i) );
end
AlphaMv(nt)=( (25-V(nt))/10 )/( exp( (25-V(nt))/10 )-1 );    BetaMv(nt)=4*exp(-V(nt)/18);
AlphaNv(nt)=(10-V(nt))/( 100*( exp( (10-V(nt))/10 )-1) );    BetaNv(nt)=0.125*exp(-V(nt)/80);
AlphaHv(nt)=0.07*exp( -V(nt)/20 );                           BetaHv(nt)=1/( exp( (30-V(nt))/10 )+1 );
%V
%-----                                Approximate solution of U                                -----%
U=zeros(1,nt);   P=zeros(1,nt);   Q=zeros(1,nt);   R=zeros(1,nt);
U(nt)=0;         P(nt)=0;         Q(nt)=0;         R(nt)=0;

for i=1:nt-1
if( V(nt-i+1)-V(nt-i)~=0 )
PsiM   =[ ( 1-m(nt-i+1) )*( AlphaMv(nt-i+1)-AlphaMv(nt-i) ) - m(nt-i+1)*( BetaMv(nt-i+1)-BetaMv(nt-i) ) ]/[ V(nt-i+1)-V(nt-i) ];
PsiN   =[ ( 1-n(nt-i+1) )*( AlphaNv(nt-i+1)-AlphaNv(nt-i) ) - n(nt-i+1)*( BetaNv(nt-i+1)-BetaNv(nt-i) ) ]/[ V(nt-i+1)-V(nt-i) ];
PsiH   =[ ( 1-h(nt-i+1) )*( AlphaHv(nt-i+1)-AlphaHv(nt-i) ) - h(nt-i+1)*( BetaHv(nt-i+1)-BetaHv(nt-i) ) ]/[ V(nt-i+1)-V(nt-i) ];
end
if( V(nt-i+1)-V(nt-i)==0 )
PsiM=0;  PsiN=0;  PsiH=0;  
end
U(nt-i)=U(nt-i+1)-d*( Gna*m(nt-i+1)^a*h(nt-i+1)^b+Gk*n(nt-i+1)^c +Gl )*U(nt-i+1) ...
                 -d*( PsiM*P(nt-i+1)+PsiN*Q(nt-i+1)+PsiH*R(nt-i+1) ) ... 
                 -d*( Vp(nt-i+1)-V(nt-i+1) );
P(nt-i)=P(nt-i+1)-( AlphaMv(nt-i+1)+BetaMv(nt-i+1) )*P(nt-i+1)*dt + a*Gna*m(nt-i+1)^(a-1)*h(nt-i+1)^b*(V(nt-i+1)-Ena)*U(nt-i+1)*dt; 
Q(nt-i)=Q(nt-i+1)-( AlphaNv(nt-i+1)+BetaNv(nt-i+1) )*Q(nt-i+1)*dt + c*Gk *n(nt-i+1)^(c-1)*(V(nt-i+1)-Ek)*U(nt-i+1)*dt; 
R(nt-i)=R(nt-i+1)-( AlphaHv(nt-i+1)+BetaHv(nt-i+1) )*R(nt-i+1)*dt + b*Gna*m(nt-i+1)^a*h(nt-i+1)^(b-1)*(V(nt-i+1)-Ena)*U(nt-i+1)*dt; 
end

