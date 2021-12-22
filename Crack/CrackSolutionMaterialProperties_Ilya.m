function [Cd, Cs, Cr, nu, ro, E, mu, Gamma, PlaneStrain, tau_p, Xc0]=CrackSolutionMaterialProperties

% % %---PMMA
Cd=2700;%should be plane strain velocity
Cs=1345;
nu=(2*(Cs/Cd).^2-1)/2/((Cs/Cd).^2-1);%plane strain
ro=1.17E3;%[Kg/m^3];
E=Cs^2*2*ro*(1+nu);%[Pa]
mu=Cs^2*ro;
%-------Choose true for plane strain in boundary conditions
PlaneStrain=false;
if(PlaneStrain~=true)
    %-------Plain Stress
    Cd=(E/ro/(1-nu^2))^0.5; %approx. =2340m/s;
end

%Gamma=10;%[J/m^2] Fast lubricated
%tau_p=2.6578e6;
Gamma=2.5;%[J/m^2] Dry
%tau_p=1.6175e6;
% 
% %%%assume Xc0
% %Xc0=4.5e-3;% for Dry
Xc0=2.7e-3;% for Dry
 %Xc0=4e-3;%  for lubricated
% 
% %Gamma=3.2;%[J/m^2] Dry
% %Xc0=2.5e-3;%
% %Gamma=2;%[J/m^2] Dry
% %Xc0=3e-3;%
% 

%----Linear cohesive zone
% 
% tau=@(x) 1+x;
% ILowerBoundery=-1;
%---
tau=@(x) exp(x);
ILowerBoundery=-10;

A=1;
k=Cs/Cd;
K=(Gamma*mu*4*(1-k^2)./A).^0.5;%Broberg p.334,336
f=@(x) tau(x).*(-x).^-0.5;
tau_p=K/(2/pi)^0.5/(integral(f,ILowerBoundery,0)*Xc0^0.5);%Eq.16c



 %---------Just for check
% Cd=2498;%should be plane strain velocity
% Cs=1200;
% nu=(2*(Cs/Cd).^2-1)/2/((Cs/Cd).^2-1);%plane strain
% ro=1.23E3;%[Kg/m^3];
% mu=Cs^2*ro;
% E=Cs^2*2*ro*(1+nu);%[Pa]
% % 
 %E=;
%  ro=1.17E3;%[Kg/m^3];
%  mu=E/2/(1+nu)*1.2;
%  Cs=(mu/ro)^0.5;
%  nu=nu*1.1;
%  E=Cs^2*2*ro*(1+nu);%[Pa]
%  Cd1=( E*(1-nu)/ro/(1+nu)/(1-2*nu) ) ^0.5;
%  Cd=(E/ro/(1-nu^2))^0.5; 

% %-------Simualtion
% 
%  E=5.65*1e9;
%  nu=0.33;
%  ro=1.18E3;%[Kg/m^3];
%  Gamma=1.12;%[J/m^2]
%  mu=E/2/(1+nu);
%  Cs=(mu/ro)^0.5;
%  Cd=( E*(1-nu)/ro/(1+nu)/(1-2*nu) ) ^0.5;
%  tau_p=1e6;%[Pa]
% % %-------Choose true for plane strain in boundary conditions
%  PlaneStrain=false;
%  if(PlaneStrain~=true)
%     %-------Plain Stress
%      Cd=(E/ro/(1-nu^2))^0.5; %approx. =2340m/s;
%  end
%  
 %------Homolite
%  Cd=2590;%should be plane strain velocity
% Cs=1284;
% nu=(2*(Cs/Cd).^2-1)/2/((Cs/Cd).^2-1);%plane strain
% ro=1.19E3;%[Kg/m^3];
% E=Cs^2*2*ro*(1+nu);%[Pa]
% mu=Cs^2*ro;
% %-------Choose true for plane strain in boundary conditions
% PlaneStrain=true;
% if(PlaneStrain~=true)
%     %-------Plain Stress
%     Cd=(E/ro/(1-nu^2))^0.5; %approx. =2340m/s;
% end
%  
 
 %-------McLasky

 %E=43.2e9;
% nu=0.112;
% ro=2670;%[Kg/m^3];
% G=0.3;%[J/m^2]
% mu=E/2/(1+nu);
% Cs=(mu/ro)^0.5;
% Cd=( E*(1-nu)/ro/(1+nu)/(1-2*nu) ) ^0.5;
% %-------Choose true for plane strain in boundary conditions
% PlaneStrain=true;
% if(PlaneStrain~=true)
%    %-------Plain Stress
%     Cd=(E/ro/(1-nu^2))^0.5; %approx. =2340m/s;
% end


%-----------------calc Cr
R=@(z) (Cs^-2 - 2*z.^2).^2 + 4*z.^2 .* (Cd^-2 -z.^2).^0.5.* (Cs^-2-z.^2).^0.5; %just after eq. 6.3.42
D=@(v) -v.^4.*real(R(1./v)); %Due to numerical error R has an imaginary part
Cr=fzero(D,Cs);


