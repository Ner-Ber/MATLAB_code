function xc=Crack_SupershearCalcXc(v,Cs,tau,ILowerBoundery)
%Calculates the size of the linear cohesive zone for super shear ruptures.  Broberg p354 6.3.70
%[v]=Cd,[Cs]=Cd
%[xc]=G*mu/tau_s^2;

tau=@(w)(1+w);
ILowerBoundery=-1;

G=1;%[J/m^2]
mu=1;% mu=Cs^2*ro;
tau_s=1;%mu/1000;
Cd=1;

k=(Cs/Cd);
b=v/Cd;%p332
ap=(1-b^2)^0.5;
bs=((b/k)^2-1)^0.5;
g=1/pi*atan(4*ap*bs/(bs^2-1)^2);
YII=(1-k^2)*b^2*sin(pi*g)/(2*k^2*(1-b^2)^0.5);


dz=1E-3;%[m]
z=dz:dz:1-dz; % the integral of I is singular at z=1;

%---------Linear cohesive zone
%6.3.70 has two internal integrals,named here by I1  and I2

% for linear cohesive zone (otherwise this hould be an integral)
%I1=1/g;

%calculate I1 for general cohesive zone %Need to check accurecy at various Cf

I1=zeros(1,length(z));

for j=1:length(z)
    I1(j)=calcI1sym(tau,ILowerBoundery,-z(j),g);
end

%calculate I2
I2=zeros(1,length(z));
for j=1:length(z)
    I2_integrand=@(w) 1./(w-z(j))./((w).^(1-g));
    I2(j)=integral(I2_integrand,1,1E7);
end

wd_integrand=tau(-z).*z.^(1-g).*(I1+tau(-z).*I2);
wd=dz*trapz(wd_integrand);

xc=G*mu/tau_s^2* (pi*(1-k^2))/(YII*sin(pi*g)*wd);%p.354 6.3.39


function I1=calcI1sym(tau,ILowerBoundery,z,g)
%Integral is done in non dimensional coordinats

I=@(w) -(tau(z)- tau(w))./(w-z)./((-w).^(1-g));

%------Decomposition to sym and Asym
v=@(w) 1/2*( I(w)+ I(2*z-w));
%h=@(w) 1/2*( I(w)- I(2*z-w));
l=min(z-ILowerBoundery,0-z);
%I1=integral(v,z-l/2,z+l/2)+integral(I,-1,z-l/2)+integral(I,z+l/2,0,'RelTol',1e-3,'AbsTol',1e-5);
%I1=integral(I,-1,z-l/2)+2*integral(v,z-l/2,z)+my_integralSingular(I,1-g,z+l/2,0);
I1=integral(I,ILowerBoundery,z-l/2)+2*integral(v,z-l/2,z)+my_integralSingular(I,1-g,z+l/2,0);

%I1=( I1-1i*pi*tau(z)/(-z)^(1-g) );%Second part in (...) needed for Principal value.


