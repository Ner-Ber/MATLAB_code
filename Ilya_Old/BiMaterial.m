function f=BiMaterial(v)

%Units:  [v]=[m/s]

%Cd=2340;

%-------Material-1
Cd1=2333;
Cs1=1345;
ro1=1.17E3;%[Kg/m^3];

nu1=(2*(Cs1/Cd1).^2-1)/2/((Cs1/Cd1).^2-1);
E1=Cs1^2*2*ro1*(1+nu1);%[Pa]
mu1=Cs1^2*ro1;

%Cd=(E/ro/(1-nu^2))^0.5;
alpha_d1=(1-(v/Cd1).^2).^0.5;
alpha_s1=(1-(v/Cs1).^2).^0.5;
D1=4*alpha_d1.*alpha_s1-(1+alpha_s1.^2).^2;

%-------Material-2
% %Cd2=(E1/ro1/(1-nu1^2))^0.5;
% Cd2=2622;
% Cs2=1473;
% ro2=1.17E3;%[Kg/m^3];
% 
% nu2=(2*(Cs2/Cd2).^2-1)/2/((Cs2/Cd2).^2-1);
% E2=Cs2^2*2*ro2*(1+nu2);%[Pa]
% mu2=Cs2^2*ro2;

nu2=0.37;
mu2=1.1*mu1;
ro2=ro1;

Cs2=(mu2/ro2)^0.5;
Cd2=(2*mu2/ro2/(1-nu2))^0.5;

alpha_d2=(1-(v/Cd2).^2).^0.5;
alpha_s2=(1-(v/Cs2).^2).^0.5;
D2=4*alpha_d2.*alpha_s2-(1+alpha_s2.^2).^2;


f=(1-alpha_s1.^2).*alpha_d1.*mu2.*D2+(1-alpha_s2.^2).*alpha_d2.*mu1.*D1;

