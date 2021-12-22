function Xc=CrackSolution_CalcXc(v,Cd)
%v= [Cs]
%Cd=[Cs]
%Xc is given in units of Gamma*mu/tau_p

Cs=1;
k=1/Cd;
Gamma=1;
mu=1;
tau_p=1;

alpha_d=(1-(v/Cd).^2).^0.5;
alpha_s=(1-(v/Cs).^2).^0.5;
D=4*alpha_d.*alpha_s-(1+alpha_s.^2).^2;
A=2*(1-k^2)*alpha_s.*v.^2./D/Cs^2;

%---------Linear slip weakening
%   tau=@(x)1-(-x).^1.4;
%   ILowerBoundery=-1;
%------
 tau=@(x) exp(x);
 ILowerBoundery=-10;

K=(Gamma*mu*4*(1-k^2)./A).^0.5;%Broberg p.334,336
f=@(x) tau(x).*(-x).^-0.5;
Xc=( K/tau_p/(2/pi)^0.5/ (integral(f,ILowerBoundery,0)) )^2;%Eq.16c