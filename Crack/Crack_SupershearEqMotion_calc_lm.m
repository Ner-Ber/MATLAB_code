function [tau_0,l_m,c_m,g_m]=Crack_SupershearEqMotion_calc_lm(c,g,Gmy,k)
%calculates the minimal length of supershear propagation

%---paramteres for k=0.48
tau_0=linspace(0.05,0.44,50); %in units of tau_d
%tau_0=0.44;
%tau_0=0.1;
[~, index]=min(abs(c-0.54));
c=c(index:end);
g=g(index:end);
Gmy=Gmy(index:end);

for j=1:length(tau_0)
a=1./(tau_0(j).^(1./g).*Gmy);
%plot(c(2:end),smooth(diff(a)./diff(c),7));
%plot(c,a);
%hold all;
[~,index]=min(abs(smooth(diff(a)./diff(c),7)));
l_m(j)=a(index);
c_m(j)=c(index);
g_m(j)=g(index);
end
