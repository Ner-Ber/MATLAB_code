function Sxy=calc_dynamic_stress(Uxy,Uxy0)
%Uxy -  strain
%Uxy0 -  initial strain

if (nargin<2)
    Uxy0=mean(Uxy(50:150,:),1);
end

Uxy=Uxy*1e-3;
Uxy0=Uxy0*1e-3;


[Cd, Cs, Cr, ~ , ~ ,E, mu, Gamma, ~ ,~]=CrackSolutionMaterialProperties;

alpha=0.43;%alpha=1-mu_inf/mu_0
%Sxy=2*mu_inf*Uxy_0+2*mu_0*(Uxy-Uxy0)

if (size(Uxy,1)==size(Uxy0,1))
%For Uxy and Uxy0 of same size
    Sxy=2*mu*(Uxy-alpha*Uxy0);
else
%For temporal matice of Uxy 
Sxy=2*mu*(Uxy-alpha*repmat(Uxy0,length(Uxy(:,1)),1) );    
end
Sxy=Sxy*1e-6;
