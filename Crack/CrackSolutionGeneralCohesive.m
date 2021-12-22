function sol=CrackSolutionGeneralCohesive(v,h)
%v[Cr],h[m]
% all formulas taken from Samudrala@Samudrala,Huang&Rosakis JGR2002

[Cd, Cs, Cr, nu , ~, ~, mu, Gamma, PlaneStrain, tau_p]=CrackSolutionMaterialProperties;

%-- manual input:
Gamma = 4;

k=Cs/Cd;%Broberg p.330

v=v*Cr;
%v=1100;
dx=0.2*10e-5;%10E-5;%[m]
%x=5E-3:-dx:-40E-3;%[m] distance from the crack tip ,theta=0;
%x=-30e-3;
%x=60E-3:-dx:-60E-3;%[m] distance from the crack tip ,theta=0;
% x=50E-3:-dx:-50E-3;%[m] distance from the crack tip ,theta=0;
x=80E-3:-dx:-80E-3;%[m] distance from the crack tip ,theta=0;
%------- General LEFM functions
alpha_d=(1-(v/Cd).^2).^0.5;
alpha_s=(1-(v/Cs).^2).^0.5;
D=4*alpha_d.*alpha_s-(1+alpha_s.^2).^2;
A=2*(1-k^2)*alpha_s.*v.^2./D/Cs^2; %Broberg p.334,336
%A=alpha_s.*v.^2./D/(1-nu)/Cs^2; %Only for plain strain

%------define complex variables
z_d=x+1i*alpha_d*h;
z_s=x+1i*alpha_s*h;

%-------Create \tau profile

%----Linear cohesive zone

%   tau=@(x) 1+x;
%   ILowerBoundery=-1;

%---------Linear slip weakening
   %tau=@(x)1-(-x).^1.4;
   %ILowerBoundery=-1;
% %------
 tau=@(x) exp(x);
 ILowerBoundery=-10;
%-------
% tau=@(x) (1+log(1-10*x)).*exp(x);
% ILowerBoundery=-10;
%-------
% tmp=load('C:\Users\owner\Documents\MATLAB\tmp\c4.mat');
% p=tmp.p;
% % tmp=load('C:\Users\owner\Documents\MATLAB\tmp\c2.mat');
% % p2=tmp.p;
% ILowerBoundery=-1;
% tau=@(x) fnval(p,x);


%---------Calc the potentials
%Note that unlike eq.13 L is not necessary ILowerBoundery.

%-------Assume Gamma and L
%  Gamma=0.3;%[J/m^2]
%  L=25E-3;

%-------Assume Gamma and L_0(Cf=0)
  L0=2.5E-3;
  L=L0/A;
% % % % 
% 
  K=(Gamma*mu*4*(1-k^2)./A).^0.5;%Broberg p.334,336
  f=@(x) tau(x).*(-x).^-0.5;
  tau_p=K/(2/pi)^0.5/(integral(f,ILowerBoundery,0)*L^0.5);%Eq.16c
%-------------------------------------------

%-------Assume Gamma and tau_p

%   Gamma=1.12;%[J/m^2]
%   tau_p=1.6e6;
%   K=(Gamma*mu*4*(1-k^2)./A).^0.5;%Broberg p.334,336
%   f=@(x) tau(x).*(-x).^-0.5;
%   L=( K/tau_p/(2/pi)^0.5/ (integral(f,ILowerBoundery,0)) )^2;%Eq.16c
%  
%--------------------------------

A0=-1i/mu/(2*pi)^0.5*2*alpha_s/D*K;%Eq.15
f=@(x,z) (-x).^0.5.*tau(x)./(x-z);

%------calc F
z=z_d;
I=tau_p*calcI(f,z/L,ILowerBoundery)*L^0.5;
F=A0*z.^-0.5+(2*pi*1i)^-1*4*alpha_s/mu/D*z.^-0.5.*I;%Eq.13

%------calc G
z=z_s;
I=tau_p*calcI(f,z/L,ILowerBoundery)*L^0.5;
F_tmp=A0*z.^-0.5+(2*pi*1i)^-1*4*alpha_s/mu/D*z.^-0.5.*I;
G=-0.5*(1+alpha_s^2)/alpha_s*F_tmp;%Eq.11

%-----Calc Stress/Strain
sol=calcStress(mu*F,mu*G,alpha_d,alpha_s);
[sol.Uxy sol.Uxx sol.Uyy]=calcStrainFromStress(sol,mu,k);

sol.x=x;
sol.y=h;
sol.v=v;
sol.t=-sol.x/sol.v;
sol.Cr=Cr;
sol.Cd=Cd;
sol.Cs=Cs;
sol.Xc=L;
sol.nu=nu;
sol.mu=mu;
sol.Gamma=Gamma;
sol.tau_p=tau_p;
sol.dc=2*sol.Gamma/sol.tau_p;
sol.smtT=ceil(1E-6/(abs(mean(diff(sol.x)))/sol.v));
sol.Ux=-cumsum(sol.Uxx)*(sol.x(1)-sol.x(2));% [m] for slip multiply by 2.

[~, index]=min(abs(sol.x*1e3+40));
sol.Uxy_Offset=sol.Uxy(index);

function I=calcI(f,z,ILowerBoundery)

for j=1:length(z)
    f_tmp=@(x) f(x,z(j));
    I(j)=integral(f_tmp,ILowerBoundery,0);
end


function o=calcStress(F,G,alpha_d,alpha_s)
o.Sxx=real((1-alpha_s^2+2*alpha_d^2)*F+2*alpha_s*G);
o.Syy=-real((1+alpha_s.^2)*F+2*alpha_s*G);
o.Sxy=-imag(2*alpha_d*F+(1+alpha_s^2)*G);

function [Uxy Uxx Uyy]=calcStrainFromStress(sol,mu,k)

%-------Works for both boundary conditions with appropriate k 

Uxy=1/(2*mu)*sol.Sxy;
Uxx=1/(2*mu)/(2*(1-k^2))*(sol.Sxx-(1-2*k^2)*sol.Syy);
Uyy=1/(2*mu)/(2*(1-k^2))*(sol.Syy-(1-2*k^2)*sol.Sxx);


