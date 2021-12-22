function sol=CrackSolutionSP(v,R,L,h)
%v[Cr],R[m],L[R0],h[m] 
% all formulas taken from Rice 2005,Appendeix
%Returns LEFM solution, and slip pulse   

Cd=2700;
Cs=1345;
nu=(2*(Cs/Cd).^2-1)/2/((Cs/Cd).^2-1);%plane stress
ro=1.17E3;%[Kg/m^2];
E=Cs^2*2*ro*(1+nu);%[Pa]
G=1;%[J/m^2]


%-------calc C_Rayleigh 
tmp=CrackSolution_n(900:1350,0,1);
[~,index]=min(abs(tmp.D));
Cr=tmp.v(index);

v=v*Cr;
dx=10E-6;%[m]
x=-70E-3:dx:70E-3;%[m] distance from the crack tip ,theta=0;

%-------
alpha_d=(1-(v/Cd).^2).^0.5;
alpha_s=(1-(v/Cs).^2).^0.5;
D=4*alpha_d.*alpha_s-(1+alpha_s.^2).^2;
A=alpha_s.*v.^2./D/(1-nu)/Cs^2; %------Note!Doesn't go to 1 when c_f->0. 
%K=(G*E/(1-nu^2)./A).^0.5; %plane strain
K=(G*E./A).^0.5; %plane stress
K0=(G*E).^0.5; %plane stress

z_d=x+1i*alpha_d*h;
z_s=x+1i*alpha_s*h;

%------LEFM
z=z_d;
M_d=K./(2*pi*z).^0.5;
z=z_s;
M_s=K./(2*pi*z).^0.5;

sol.Lefm=calcStress(M_d,M_s,alpha_d,alpha_s,D);
[sol.Lefm.Uxy sol.Lefm.Uxx sol.Lefm.Uyy sol.Lefm.UxxPstrn sol.Lefm.UyyPstrn]=calcStrainFromStress(sol.Lefm,nu,E);

%-----Slip Pulse

R0=2.5E-3;%[m]
tau_p=3/4*(pi/2)^0.5/R0^0.5*K0; % (A12a)

R=9*pi/32*K^2/tau_p^2;
%R=R0; %Good only for slow fronts
L=L*R;

Theta_tag=2*asin((R/L)^0.5);
M0=-1/pi*tau_p*( Theta_tag-(Theta_tag-sin(Theta_tag))/(2*sin(Theta_tag/2)^2) );

F1=( exp(1i*Theta_tag)-1 )*L/2;
F2=@(z) z.^0.5.*(z+L).^0.5 ;
F=@(z) (F1-z-F2(z))./(F1-z+F2(z)).*(z-F2(z))./(z+F2(z));

M1=@(z) -1/pi*tau_p*( (1+z/R)*(-1i).*log(F(z))+z.^0.5.*(z+L).^0.5*Theta_tag/R);


sol.Sp=calcStress(M0+M1(z_d),M0+M1(z_s),alpha_d,alpha_s,D);
[sol.Sp.Uxy sol.Sp.Uxx sol.Sp.Uyy sol.Sp.UxxPstrn sol.Sp.UyyPstrn]=calcStrainFromStress(sol.Sp,nu,E);
sol.Sp.Xc=R;
sol.Sp.L=L;
sol.Sp.tau_p=tau_p;
sol.x=x;
sol.v=v;
sol.Cr=Cr;

function O=calcStress(M_d,M_s,alpha_d,alpha_s,D)

O.Sxx=2*alpha_s./D.*imag((1+2*alpha_d.^2-alpha_s.^2).*M_d-(1+alpha_s.^2).*M_s);
O.Syy=-2*alpha_s.*(1+alpha_s.^2)./D.*imag(M_d-M_s);
O.Sxy=1./D.*real(4*alpha_s.*alpha_d.*M_d-(1+alpha_s.^2).^2.*M_s);

function [Uxy Uxx Uyy UxxPstrn UyyPstrn]=calcStrainFromStress(sol,nu,E)
%planes stress
Uxy=(1+nu)/E*sol.Sxy;
Uxx=1/E*(sol.Sxx-nu*sol.Syy);
Uyy=1/E*(sol.Syy-nu*sol.Sxx);
%planes strain
UxxPstrn=1/E*(1+nu)*(sol.Sxx*(1-nu)-nu*sol.Syy);
UyyPstrn=1/E*(1+nu)*(sol.Syy*(1-nu)-nu*sol.Sxx);

