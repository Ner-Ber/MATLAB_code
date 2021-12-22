function [Gamma]=CalcGamma%(L,R)

% vg=60/(R.x0-L.x0);
% 
% Gamma.vg=vg;
% Gamma.w=R.w;
% 
% %Calc From Sigma
% Gamma.Gamma(1)=(R.s^2-L.s^2)^0.5*L.s/(R.x0-L.x0)*vg^2;
% 
% %Calc From Amplitude
% 
% %Gamma.Gamma(2)=((L.a/1.3/R.a)^2-1)^0.5*(L.s*vg)^2/(R.x0-L.x0);
% Gamma.Gamma(2)=((L.a/1.3/R.a)^4-1)^0.5*(L.s*vg)^2/(R.x0-L.x0);
% %theory
% 
% Gamma.Gamma(3)=L.w/4;

%-stam

w=6e6;
wc=10e6;
beta=2*asin(w/wc);
vg=wc/2*cos(beta/2);
t=60/vg;
G=w/4;
s0t=2.2e-6;
s0=s0t*vg;
s=(s0^2+(G*t/s0)^2)^0.5;
A=(s0/s)^0.5;

%---
s0t=2.2e-6;
s0w=1/2.2e-6;
Dbeta=beta-2*asin((w+0.5*s0w)/wc);
s0x=1/2/Dbeta;

%------
a4=1/(2*3)*wc/8*cos(beta/2);

coeffRatio=a4*Dbeta/(1/2*G);
A