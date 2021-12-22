function sol=CrackSolutionSW(v,h)
%v[Cr],h[m]
% all formulas taken from Poliakov&Rice 2002,Appendeix

[Cd Cs Cr nu ro E mu G PlaneStrain]=CrackSolutionMaterialProperties;
k=Cs/Cd;%Broberg p.330

v=v*Cr;
dx=10E-6;%[m]
x=-70E-3:dx:70E-3;%[m] distance from the crack tip ,theta=0;

%-------
alpha_d=(1-(v/Cd).^2).^0.5;
alpha_s=(1-(v/Cs).^2).^0.5;
D=4*alpha_d.*alpha_s-(1+alpha_s.^2).^2;
A=2*(1-k^2)*alpha_s.*v.^2./D/Cs^2; %Broberg p.334,336

z_d=x+1i*alpha_d*h;
z_s=x+1i*alpha_s*h;

%--Assume G and tau_p
% tau_p=3.16e5;%maximal shear stress
% G=0.3;
% K=(G*mu*4*(1-k^2)./A).^0.5;%Broberg p.334,336
% R=9*pi/32*K.^2/tau_p^2;

%--Assume G and R
R=0.025;%[m] Cohesive zone size
G=0.3;
K=(G*mu*4*(1-k^2)./A).^0.5;%Broberg p.334,336
tau_p=(9*pi/32*K.^2/R)^0.5;

%---------
M=@(z) 2/pi.*tau_p.*( (1+z/R).*atan((z/R).^-0.5)-(z/R).^0.5 ); %(A11)

sol=calcStress(M(z_d),M(z_s),alpha_d,alpha_s,D);
[sol.Uxy sol.Uxx sol.Uyy]=calcStrainFromStress(sol,mu,k);
sol.Xc=R;
sol.x=x;
sol.v=v;
sol.Cr=Cr;
sol.t=-sol.x/sol.v;

function O=calcStress(M_d,M_s,alpha_d,alpha_s,D)

O.Sxx=2*alpha_s./D.*imag((1+2*alpha_d.^2-alpha_s.^2).*M_d-(1+alpha_s.^2).*M_s);
O.Syy=-2*alpha_s.*(1+alpha_s.^2)./D.*imag(M_d-M_s);
O.Sxy=1./D.*real(4*alpha_s.*alpha_d.*M_d-(1+alpha_s.^2).^2.*M_s);

function [Uxy Uxx Uyy]=calcStrainFromStress(sol,mu,k)

%-------Works for both boundary conditions with appropriate k 

Uxy=1/(2*mu)*sol.Sxy;
Uxx=1/(2*mu)/(2*(1-k^2))*(sol.Sxx-(1-2*k^2)*sol.Syy);
Uyy=1/(2*mu)/(2*(1-k^2))*(sol.Syy-(1-2*k^2)*sol.Sxx);

