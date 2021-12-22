function sol = CrackSolutionGeneralCohesive_InsertVariables(v,h,Gamma,variableName1,variableVal1,varargin)
% sol = CrackSolutionGeneralCohesive_InsertVariables(v,h,Gamma,variableName1,variableVal1,cohesiveFun,PlaneFlag,x_range,dx)
%
% CrackSolutionGeneralCohesive_InsertVariables will calculate the cohesive
% zone model in for velocity 'v' at hight 'h'.
%
% INPUTS:
%v - velocity in units of [Cr]
%h - height in units of [m]
%Gamma - fracture energy [J/m^2]
% variableName1,variableVal1 - is the variable by which the cohesive model
% is calculated, for example the slip length 'L0' and its value '25e-3'.
% The optional variables are 'L0','L','tau_p'.
%
% OPTIONAL:
% cohesiveFun - function describing the cutoff of the cohesive zone:
%               options: 'exp' (default), 'linear', 'Ohnaka89', 'slipWeak'
% PlaneFlag - define as 'PlaneStrain' (default) or 'PlaneStress'
% x_range - spacial vector definning solution range (defaults: 80E-3:-dx:-80E-3 [m])
% dx - spacial vector x spacing. (defaults: 5e-6 [m])
%
% all formulas taken from Samudrala@Samudrala,Huang&Rosakis JGR2002

%% set defaults and parameters
dx_def = 5e-6;   %[m]
x_range_def = [80E-3,-80E-3]; %[m] distance from the crack tip ,theta=0;
[cohesiveFun,PlaneFlag,x_range,dx] = setDefaults4function(varargin,'exp','PlaneStrain',x_range_def,dx_def);

x = max(x_range):-abs(dx):min(x_range); %[m] distance from the crack tip ,theta=0;

[Cd, Cs, Cr, nu , ~, ~, mu, GammaDefault, PlaneStrain, tau_p_def]=CrackSolutionMaterialProperties(PlaneFlag);

k=Cs/Cd;%Broberg p.330
v=v*Cr;

%------- General LEFM functions
alpha_d=(1-(v/Cd).^2).^0.5;
alpha_s=(1-(v/Cs).^2).^0.5;
D=4*alpha_d.*alpha_s-(1+alpha_s.^2).^2;

if ~strcmpi(PlaneFlag,'PlaneStrain')
    A=2*(1-k^2)*alpha_s.*v.^2./D/Cs^2; %Broberg p.334,336, equation 6.2.26
else
    A=alpha_s.*v.^2./D/(1-nu)/Cs^2; %Only for plain strain
end

%------define complex variables
z_d=x+1i*alpha_d*h;
z_s=x+1i*alpha_s*h;

%% Create \tau profile

if strcmpi(cohesiveFun,'linear') %-- Linear cohesive zone
    tau=@(x) 1+x;
    ILowerBoundery=-1;
    
elseif strcmpi(cohesiveFun,'slipWeak') %-- Linear slip weakening
    tau=@(x)1-(-x).^1.4;
    ILowerBoundery=-1;
    
elseif strcmpi(cohesiveFun,'exp') %-- exponential cohesive
    tau=@(x) exp(x);
    ILowerBoundery=-10;
    
elseif strcmpi(cohesiveFun,'Ohnaka89') %---- Cohesive model from Ohnaka&Yamashita JGR 89
    tau=@(x) (1+log(1-10*x)).*exp(x);
    ILowerBoundery=-10;
    
elseif isa(cohesiveFun,'function_handle')
    tau= cohesiveFun;
%     ILowerBoundery=-10;
    ILowerBoundery=-1;
    
elseif loadCohesive %-------
    % tmp=load('C:\Users\owner\Documents\MATLAB\tmp\c4.mat');
    % p=tmp.p;
    % % tmp=load('C:\Users\owner\Documents\MATLAB\tmp\c2.mat');
    % % p2=tmp.p;
    % ILowerBoundery=-1;
    % tau=@(x) fnval(p,x);
end


%% Calc the potentials
%Note that unlike eq.13 L is not necessary ILowerBoundery.

if strcmpi(variableName1,'L')   %-- Assume Gamma and L
    %     L=25E-3;
    L = variableVal1;
    K=(Gamma*mu*4*(1-k^2)./A).^0.5;%Broberg p.334,336 Eq. 6.2.42
    tauDecay=@(x) tau(x).*(-x).^-0.5;
    tau_p=K/(2/pi)^0.5/(integral(tauDecay,ILowerBoundery,0)*L^0.5);%Eq.16c
    
elseif strcmpi(variableName1,'L0')   %-- Assume Gamma and L0 (Cf=0)
    %     L0=2.5E-3;
    L0 = variableVal1;
    L=L0/A;                         % is this equation 46? (Neri)
    K=(Gamma*mu*4*(1-k^2)./A).^0.5;%Broberg p.334,336 Eq. 6.2.42
    tauDecay=@(x) tau(x).*(-x).^-0.5;
    tau_p=K/(2/pi)^0.5/(integral(tauDecay,ILowerBoundery,0)*L^0.5);%Eq.16c
    
elseif strcmpi(variableName1,'tau_p')   %-- Assume Gamma and tau_p
    %     tau_p=1.6e6;
    tau_p = variableVal1;
    K=(Gamma*mu*4*(1-k^2)./A).^0.5;%Broberg p.334,336 Eq. 6.2.42
    tauDecay=@(x) tau(x).*(-x).^-0.5;
    L=( K/tau_p/(2/pi)^0.5/ (integral(tauDecay,ILowerBoundery,0)) )^2;%Eq.16c
end
%% calc solution

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
sol.tauDecayOfx = tauDecay;
sol.tauOfx = tau;

[~, index]=min(abs(sol.x*1e3+40));
sol.Uxy_Offset=sol.Uxy(index);

%% external functions

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


