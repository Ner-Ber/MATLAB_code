function [Uxy_0,vO]=Crack_SupershearEqMotionCalcV2(sol,tau_0)
%Changing profile of the interface paramters
% use [sol]=Crack_SupershearEqMotion2 as an imput

[~, Cs ,~, ~ ,~ ,~, mu, Gamma ,~,tau_d, Xc0]=CrackSolutionMaterialProperties;

Cd=Cs/sol.k;



%For x<100mm
Gamma1=Gamma/1.65;
tau_d1=tau_d/1.65^0.5;
a01=Gamma1/tau_d1^2*mu; %Xc doesn't change
tau_d1=1.3e6;
%For x>100
Gamma2=Gamma;
tau_d2=tau_d;
a02=Gamma2/tau_d2^2*mu;
tau_d2=1.6e6;

%For x>130mm
% Gamma3=Gamma;
% tau_d3=tau_d;
% a03=Gamma3/tau_d3^2*mu;

%tau_0=linspace(0.05*tau_d,0.98*tau_d,20);
%Uxy_0=tau_0/2/mu;

%Uxy_0=0.225e-3;
%Uxy_0=[0.14 0.18 0.24]*1e-3;
%tau_0=2*mu*Uxy_0;


for j=1:length(tau_0)
    a1=a01./(sol.G_my.*(tau_0(j)/tau_d1).^(1./sol.g) ) ;
    v1=sol.b*Cd;
    a1=a1(v1>1.08*Cs);
    v1=v1(v1>1.08*Cs);
    [~,index]=min(a1);
    a1=a1(index:end);
    v1=v1(index:end);
    atmp=(00:0.1:100)*1e-3;
    v1=interp1(a1,v1,atmp);
    a1=atmp;
    
    
    a2=a02./(sol.G_my.*(tau_0(j)/tau_d2).^(1./sol.g) ) ;
    v2=sol.b*Cd;
    a2=a2(v2>1.08*Cs);
    v2=v2(v2>1.08*Cs);
    [~,index]=min(a2);
    a2=a2(index:end);
    v2=v2(index:end);
    atmp=(100:0.1:200)*1e-3;
    v2=interp1(a2,v2,atmp);
    a2=atmp;
    
    %     a3=a03./(sol.G_my.*(tau_0(j)/tau_d3).^(1./sol.g) ) ;
    %     v3=sol.b*Cd;
    %     a3=a3(v3>2^0.5*Cs);
    %     v3=v3(v3>2^0.5*Cs);
    %
    %     atmp=(135:1:200)*1e-3;
    %     v3=interp1(a3,v3,atmp);
    %     a3=atmp;
    
    a=[a1  a2];
    v=[v1  v2];
    
    %-----plot
    plot(a*1e3,v,'black-');
    %my_legend_add(['k=' num2str(sol.k) ', G=' num2str(Gamma2) ', \tau_d=' num2str(round(tau_d2*1e-5)/10) ', \tau_0/\tau_p=' num2str(round(tau_0(j)/tau_d*100)/100)]);
    %my_legend_add(['k=' num2str(sol.k) ', G=' num2str(Gamma) ', \tau_0/\tau_p=' num2str(round(tau_0(j)/tau_d*100)/100), '\tau_0=' num2str(round(tau_0*1e-5)/10)]);
    my_legend_add(['k=' num2str(sol.k) ', G=' num2str(Gamma) ', Xc0=' num2str(Xc0*1e3) ', \tau_d=' num2str(round(tau_d2*1e-5)/10)  ', \tau_0=' num2str(round(tau_0*1e-4)/100) ]);
    hold all;
end