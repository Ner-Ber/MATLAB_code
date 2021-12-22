E=3;

f=0.1:10^6;
tau=10^-3;
w=2*pi*f*tau;

Eeff=E*(1+1./(1-1i./w));
close all;
semilogx(f,real(Eeff),'.-');
figure;
semilogx(f,tan(angle(Eeff)),'.-');