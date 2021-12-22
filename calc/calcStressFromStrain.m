function [Syy,Sxx,Sxy]=calcStressFromStrain(acq,PlaneStrain)

ro=1.17E3;%[Kg/m^2];
Cd=2700;
Cs=1345;
sigma=(2*(Cs/Cd).^2-1)/2/((Cs/Cd).^2-1);%plane strain

E=Cs^2*2*ro*(1+sigma);%[Pa]
% sigma=1/3;%poisson ratio
% E=3E+9/10^6 ;%youngs modulus [MPa]
% E=E/1000; %because U is in mStrain

%E=1; %used for fit

Uxx=acq.Uxx;
Uyy=acq.Uyy;
Uxy=acq.Uxy;

if(PlaneStrain == true)
    Sxx=-E/(1+sigma)/(1-2*sigma)*((1-sigma)*Uxx+sigma*Uyy);
    Syy=-E/(1+sigma)/(1-2*sigma)*((1-sigma)*Uyy+sigma*Uxx);
else
    Sxx=-E/(1-sigma^2)*(Uxx+sigma*Uyy);  % "-" is added positive stress is compresion
    Syy=-E/(1-sigma^2)*(Uyy+sigma*Uxx); % "-" is added positive stress is compresion
end

%Plane Stress (Szz=0)
% Sxx=-E/(1-sigma^2)*(Uxx+sigma*Uyy);  % "-" is added positive stress is compresion
% Syy=-E/(1-sigma^2)*(Uyy+sigma*Uxx); % "-" is added positive stress is compresion
% Sxy=E/(1+sigma)*Uxy;

%Plane Strain (Uzz=0)
% Sxx=-E/(1+sigma)/(1-2*sigma)*((1-sigma)*Uxx+sigma*Uyy);
% Syy=-E/(1+sigma)/(1-2*sigma)*((1-sigma)*Uyy+sigma*Uxx);
Sxy=E/(1+sigma)*Uxy;