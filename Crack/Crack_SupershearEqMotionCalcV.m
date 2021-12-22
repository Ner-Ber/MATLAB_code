function [Uxy_0,vO]=Crack_SupershearEqMotionCalcV(sol)
% use [sol]=Crack_SupershearEqMotion3 as an imput 

[Cd, Cs ,~, ~ ,~ ,~, mu, Gamma ,~,tau_d, Xc0]=CrackSolutionMaterialProperties;

Cd=Cs/sol.k;
a0=Gamma/tau_d^2*mu;
%tau_0=linspace(0.05*tau_d,0.98*tau_d,20);
%Uxy_0=tau_0/2/mu;

%Uxy_0=0.14e-3;
%Uxy_0=0.25e-3;
%tau_0=2*mu*Uxy_0;
tau_0=1.3e6;
%tau_0=(0.1:0.005:0.99)*tau_d;

Uxy_d=tau_d/2/mu;

for j=1:length(tau_0)
    a=a0./(sol.G_my.*(tau_0(j)/tau_d).^(1./sol.g) ) ;
    v=sol.b*Cd;

    v_cut=v(v>1.08*Cs);
    a_cut=a(v>1.08*Cs);
    [a_min(j),index]=min(a_cut);
    a_cut=a_cut(index:end);
    v_cut=v_cut(index:end);
    
    %plot(a_cut,v_cut,'o-')
    l=120e-3;
    vO(j)=interp1(a_cut,v_cut,l);
    
    %-----plot
    %plot(a*1e3,v,'black.-');
    %plot(a,(sol.b-sol.k)/(1-sol.k),'g.-');
    %my_legend_add(['k=' num2str(sol.k) ', \tau_0/\tau_p=' num2str(tau_0(j)/tau_d)]);
    %hold all;
end

plot(tau_0*1e-6,vO,'black-');
my_legend_add(['k=' num2str(sol.k) ', G=' num2str(Gamma) ', Xc0=' num2str(Xc0*1e3) ', l=' num2str(l*1e3) ]);
