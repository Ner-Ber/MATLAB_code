function [mu]=dynamicModulusPMMA(f)
%The function extrapolates the measuremnts by Read and Duncan 1981 to any
%frequnecy

%------Values are given by Read and Duncan 1981
Rmu_V=[2.28 2.29 2.24 2.14 2.09 2.08 2.07 2.05 2.02 1.99 1.96 1.95 1.86 1.76 1.66 1.54 1.52 1.5 1.47 1.42 1.38 1.32 1.26 1.21 ]*1e9;%real part of mu 
Rmu_f_V=[3597493 2529298 500035 9886 3656 2831 2113 1489 951 489 372 277 98.6 31.2 9.6 3.24 2.51 1.77 1.27 0.735 0.385 0.104 0.0329 0.0098 ];
tand_V=[0.012 0.013 0.018 0.03 0.035 0.036 0.037 0.039 0.041 0.047 0.049 0.051 0.061 0.073 0.081 0.085 0.084 0.085 0.083 0.083 0.085 0.082 0.081 0.08 0.079 0.071 ]; %imag{mu}/real{mu}
tand_f_V=[3572728 2648500 473151 9616 3388 2917 2163 1524 891.3 489 349 275 102 32.7 9.75 3.19 2.39 1.742 1.36 1.104 0.98 0.74 0.55 0.39 0.32 0.105 ];

%----sort in ascending order
[Rmu_f_V,index]=sort(Rmu_f_V);
Rmu_V=Rmu_V(index);
[tand_f_V,index]=sort(tand_f_V);
tand_V=tand_V(index);

%-------Renormalize
mu_0=2.1e9;%
mu_inf=1.2e9;
Rmu_V=(Rmu_V-Rmu_V(1))/(Rmu_V(end)-Rmu_V(1))*(mu_0-mu_inf)+mu_inf;

%------Interpolat and creat complex mu

% Rmu=@(f) interp1(log10(Rmu_f_V),Rmu_V,log10(abs(f)),'cubic');  
 %tand=@(f) interp1(log10(tand_f_V),tand_V,log10(abs(f)),'cubic'); 

Rmu=@(f) csaps(log10(Rmu_f_V),Rmu_V,0.9,log10(abs(f)));  %Note there is a problem with f=0
tand=@(f) csaps(log10(tand_f_V),tand_V,0.99,log10(abs(f))); 
Imu=@(f) Rmu(f).*tand(f);
mu=@(f) (heaviside(f).*(Rmu(f)+1j*Imu(f)) + heaviside(-f).*(Rmu(f)-1j*Imu(f)));%the real part should be simetric imaginary anti-simetric.  

mu=mu(f).';

%--To fix the problem at f=0;
 [~,index]=min(abs(f));
 mu(index)=Rmu_V(1);

%-----------Linear standard model
% mu_0=2.1e9;%mu for strain calculation
% mu_inf=1.2e9;
% tau=1/(2*pi*1500);%
% mu=@(w) (mu_inf+1j*w*tau*mu_0)./(1+1j*w*tau); %Sxy=2*mu*Uxy
% mu=mu(2*pi*f).';
