G=1;
E=5.65e9;
nu=0.33;
Kc=(G*E/(1-nu^2)).^0.5; 

%tau=0.63e6;
tau=0.2e6;
h=30e-3;
s=0.4:0.01:0.9;
a=s*h./(1-s);

Kstrip=tau*(3*h)^0.5*(1-s).^-1.*(1.13*s+(1-s)*0.28);%Tada p.299
Kstrip=Kstrip/Kc;

Kinf=tau*(pi*a).^0.5;%p.124;
Kinf=Kinf/Kc;

Kinf=1.12*tau*(pi*a).^0.5;%p.193;
Kinf=Kinf/Kc;