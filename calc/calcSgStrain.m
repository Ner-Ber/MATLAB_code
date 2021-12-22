function [U1,U2,U3]=calcSgStrain(Uxx,Uxy,Uyy)

U1=(Uxx+Uyy+2*Uxy)/2; %two bottom legs of the rossete
U3=(Uxx+Uyy-2*Uxy)/2;
U2=Uyy;

% Inverse transformation
%  Uxy=(U1-U3)/2;
%  Uyy=U2;                     %calculated normal strain
%  Uxx=U1+U3-Uyy; %calculated transvers strain