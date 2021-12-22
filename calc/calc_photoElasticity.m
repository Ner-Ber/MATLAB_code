function Iout=calc_photoElasticity(acq,theta,phi)
%Ein=[cos(theta/2);sin(theta/2)*exp(1i*phi)];

%define Parameters
C=60E-13;%[mks]
L=6E-3;
Lambda=670E-9;
ai=20/180*pi;%angle of incident
ni=1.5;
no=1;

Ein=[cos(theta/2);sin(theta/2)*exp(1i*phi)];
I=fresnelEquations(ni,no,ai);
[s1 s2 alpha]=calcPrincipalStresses(acq);% [s]=Mpa, [alpha]=deg;
alpha=alpha*pi/180;
Retardance=2*pi/Lambda*L*C*(s1-s2)*1E6; %1E6 -> convert from MPa

for j=1:length(s1(1,:)) %apply separatly for each sg
Eout=multiprod(Rot(alpha(:,j)),Ein);
Eout=multiprod(Ret(Retardance(:,j)),Eout);
Eout=squeeze(multiprod(Rot(-alpha(:,j)),Eout));


Iout.ITParl(:,j)=abs(Eout(1,:)).^2;
Iout.ITPerp(:,j)=1/2*abs(Eout(2,:)).^2;

%Transmition by Fresnel equation
% Iout.ITParl(:,j)=I.TParl.*abs(Eout(1,:)).^2;
% Iout.ITPerp(:,j)=I.TPerp.*abs(Eout(2,:)).^2;

%  Iout.ITParl(:,j)=I.TParl.*Iout.ITParl(:,j);
%  Iout.ITPerp(:,j)=I.TPerp.*Iout.ITPerp(:,j);

Iout.IT(:,j)=Iout.ITParl(:,j)+Iout.ITPerp(:,j);

%Reflection by Fresnel equations
% Iout.IRParl(:,j)=I.RParl.*abs(Eout(1,:)).^2;
% Iout.IRPerp(:,j)=I.RPerp.*abs(Eout(2,:)).^2;
% Iout.IR(:,j)=Iout.IRParl(:,j)+Iout.IRPerp(:,j);
end
Iout.t=acq.t;
 
function R=Rot(alpha)
R(1,1,:)=cos(alpha);
R(1,2,:)=-sin(alpha);
R(2,1,:)=sin(alpha);
R(2,2,:)=cos(alpha);

function R=Ret(Retardence)
R(1,1,:)=exp(1i*Retardence);
R(1,2,:)=0;
R(2,1,:)=0;
R(2,2,:)=1;