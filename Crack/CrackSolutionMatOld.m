function O=CrackSolutionMat2(v)
%v should be fraction of cr. i.e v=0.9 

E=1;
nu=1/3;

%-------calc Cr 
tmp=CrackSolution(900:1370,0);
[~,index]=min(abs(tmp.D));
Cr=tmp.v(index);

h_vec=[linspace(-7,-0.1,100)*1E-3, linspace(0.1,7,100)*1E-3];%[m]
v=v*Cr;
dx=5E-5;%[m]
x=70E-3:-dx:-70E-3;%[m]

O.x=x;
O.h=h_vec;
for j=1:length(h_vec)

h=h_vec(j);
theta=atan(h./x);

if (h>0)
%theta should be 0<theta<pi;
index=find(theta<0);
theta(index)=theta(index)+pi;
else
%theta should be -pi<theta<0;
index=find(theta>0);
theta(index)=theta(index)-pi;  
end

%smooth parameters
smtT=1E-6*v;
smtT=ceil(smtT/dx);

dThetaMax=max(diff(theta));
smtTheta=2*atan(0.5/3);
smtTheta=ceil(smtTheta/dThetaMax);

r=(h.^2+x.^2).^0.5;

sol=CrackSolution(v,theta);

O.Sxy(j,:)=sol.K.*sol.Sigma12./(2*pi*r).^0.5;
O.Sxx(j,:)=sol.K.*sol.Sigma11./(2*pi*r).^0.5;
O.Syy(j,:)=sol.K.*sol.Sigma22./(2*pi*r).^0.5;


O.Uxy=(1+nu)/E*O.Sxy;
%planes stress
O.Uxx=1/E*(O.Sxx-nu*O.Syy);
O.Uyy=1/E*(O.Syy-nu*O.Sxx);


end

figure;
% subplot(1,2,1);mesh(O.x,O.h,-O.Uxx);view(0,90);colorbar;
% subplot(1,2,2);mesh(O.x,O.h,O.Uxy);view(0,90);colorbar;
mesh(O.x*1000,O.h*1000,0.5*O.Uxy*1E-3);view(0,90);colorbar;
h = colorbar;
title(h,'U_{xy}')

