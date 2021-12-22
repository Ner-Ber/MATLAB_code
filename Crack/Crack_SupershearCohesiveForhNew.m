function sol=Crack_SupershearCohesiveForhNew(v,y)
%The function calculates elastic fields of super shear rupture according to
%Broberg p348. 
%Inseatd of specifying g_in enter veolcity [m/s].
%[y]=m

%----------Calculations done only on the lower half plane,y<0.
%If y>0 , calculations done for y<0 and then Sxx and Syy multiplied by -1.

if (y>0)
    LowerHalfPlane=false;
    y=-y; %to make sure that y <0
else
    LowerHalfPlane=true;
end
%----------------
clear global mu g ap bs k
global mu g ap bs k 

[Cd, Cs ,~, ~ ,~ ,~, mu Gamma ,~,tau_p]=CrackSolutionMaterialProperties;

%--------calc v and g
% %--Cf>2^0.5Cs
% f=@(Cf) 4./( (Cf/Cs).^2-2 ).^2 .*( ( 1-(Cf/Cd).^2 ).*(  (Cf/Cs).^2-1) ).^0.5;
% v_tmp=2^0.5*Cs:1E-4*Cs:Cd;
% [~,index]=min(abs(tan(pi*g)-f(v_tmp)));
% v=v_tmp(index);
%---------------------

k=(Cs/Cd);%general

b=v/Cd;%p332
ap=(1-b.^2).^0.5;
bs=((b/k).^2-1).^0.5;
alpha=bs/ap;
y=ap*y; %normlized units eq.6.3.22
phi_f=asin(k/b);
g=1/pi*atan(4*ap*bs/(bs^2-1)^2);
%------Frictional properties
%Unlike Sub-Rayleigh two frictional properties should be specified.

%tau=@(w) 1; %Dugdale
%tau=@(w) (1+w); %linear weakening cohesive zone in Normlized Units(w=z/Xc)
%tau=@(w) 1-(-w).^1.4;
%ILowerBoundery=-1;%Lower boundary of the integration. cohesive zone size

 tau= @(w) exp(w); %exponential 
 ILowerBoundery=-10;%Lower boundary of the integration. cohesive zone size


%---try
% Xc=5e-3; %[J/m^2]
% tau_p=1E6;
% Xc0=1;
% Gamma=1;

%------Define Gamma (fracture energy) and Xc
% 
% Gamma=3.4427;
% Xc=2.75E-3;
% Xc_n=Crack_SupershearCalcXc(v/Cd,Cs/Cd);%eq6.3.69
% tau_p=(Xc_n/Xc*Gamma*mu)^0.5;%peak shear stress
%  Xc_n0=Crack_SupershearCalcXc(2^0.5*Cs/Cd,Cs/Cd);
%  Xc0=Xc_n0*Gamma*mu/tau_p^2;

%------Define Gamma and tau_p
 %Gamma=1.5; %[J/m^2]
 %tau_p=1.4E6;%maximal shear stress
 %tau_p=2.5E6;%maximal shear stress
%
Xc_n=Crack_SupershearCalcXc(v/Cd,Cs/Cd,tau,ILowerBoundery);
Xc=Gamma*mu/tau_p^2*Xc_n;%Should be in [m]
%Xc=1.8e-3;
Xc_n0=Crack_SupershearCalcXc(2^0.5*Cs/Cd,Cs/Cd,tau,ILowerBoundery);
Xc0=Xc_n0*Gamma*mu/tau_p^2;

%-----Define Gamma and Xc at2^0.5Cs
% Gamma=1.1; %[J/m^2]
%  Xc0=7.9E-3;
%  Xc_n0=Crack_SupershearCalcXc(2^0.5*Cs/Cd,Cs/Cd);
%  tau_p=(Xc_n0/Xc0*Gamma*mu)^0.5;
%  Xc_n=Crack_SupershearCalcXc(v/Cd,Cs/Cd);
%  Xc=Gamma*mu/tau_p^2*Xc_n;%Should be in [m]

%----x and y should be normlized by Xc
dx=0.1;
%dx=0.025; %good resolution needed for the mask
%dx=0.01;

xStart=-50E-3;%[m]
xEnd=50E-3;%[m]
x=xStart:dx*Xc:xEnd;%[Xc] distance from the crack tip ,theta=0;
x=x/Xc;
y=y/Xc;


%------------Cohesive zone calculation
F=zeros(1,length(x));
f=zeros(1,length(x));
for j=1:length(x)
    
    %--------calc F_tmp for later calculation of f
    if ((x(j)-alpha*y<0 ))
        if (x(j)-alpha*y>=ILowerBoundery)
            %(x(j)-alpha*y)  %to follow the calculation
            Ftmp=calcFsym(tau,ILowerBoundery,x(j)-alpha*y);
            f(j)=-(bs^2-1)/(4*bs)*(Ftmp+conj(Ftmp));%eq.6.3.37
        else
            Ftmp=calcFNumerical(tau,ILowerBoundery,x(j)-alpha*y);
            f(j)=-(bs^2-1)/(4*bs)*(Ftmp+conj(Ftmp));
        end
    else
        f(j)=0;
        Ftmp=0;
    end
    
    %--------calc F eq.6.3.49
    %Symbolic calculation is long. When it is not necessary numerical
    %integration is used.
    if (y~=0)
        F(j)=calcFNumerical(tau,ILowerBoundery,x(j)+1i*y);
    elseif (x(j)>0)
        F(j)=calcFNumerical(tau,ILowerBoundery,x(j)+1i*y);
    else
        F(j)=Ftmp;
    end
    
end

f=f*tau_p;
F=F*tau_p;
%F=F*tau_p*0;% for Mach cone only


[soltmp.Sxx,soltmp.Syy,soltmp.Sxy]=calcStress(F,f);
%[soltmp.Sxx,soltmp.Syy,soltmp.Sxy]=calcStress(0,f);
soltmp.x=x*Xc;
soltmp.y=y/ap*Xc;

if(LowerHalfPlane==false)
    soltmp.Sxx=-soltmp.Sxx;
    soltmp.Syy=-soltmp.Syy;
    soltmp.y=y/ap*Xc;
    y=-y;
end

%---------interpolate the data at to original x
%comment those lines to avoid interp

% sol.x=xStart:dx*Xc:xEnd;
% sol.y=y/ap*Xc;
% sol.Sxx=spline(soltmp.x,soltmp.Sxx,sol.x);
% sol.Syy=spline(soltmp.x,soltmp.Syy,sol.x);
% sol.Sxy=spline(soltmp.x,soltmp.Sxy,sol.x);

sol=soltmp;

%-------Works for both boundary conditions with appropriate k 
 sol.Uxy=1/(2*mu)*sol.Sxy;
 sol.Uxx=1/(2*mu)/(2*(1-k^2))*(sol.Sxx-(1-2*k^2)*sol.Syy);
 sol.Uyy=1/(2*mu)/(2*(1-k^2))*(sol.Syy-(1-2*k^2)*sol.Sxx);
 %[sol.Uxx,sol.Uyy,sol.Uxy]=calcStrain(sol.Sxx,sol.Syy,sol.Sxy);

%------------
sol.phi_f=phi_f;
sol.g=g;
%sol.x=x*Xc;
%sol.y=y/ap*Xc;
sol.v=v;
sol.t=-sol.x/sol.v;
sol.Cs=Cs;
sol.Cd=Cd;
sol.Cr=NaN;
sol.Xc=Xc;
sol.Xc0=Xc0;
sol.Gamma=Gamma;
sol.tau_p=tau_p;
sol.smtT=ceil(1E-6/ (abs( mean( diff(sol.x))) / sol.v) );
sol.dc=2*sol.Gamma/sol.tau_p; %Approximat

Ux=cumsum(sol.Uxx(end:-1:1))*(sol.x(1)-sol.x(2)); % [m] for slip multiply by 2.
sol.Ux=Ux(end:-1:1);
sol.title=['k=' num2str(round(sol.Cs/sol.Cd*1e3)/1e3) 'Cf=' num2str(sol.v) 'G=' num2str(sol.Gamma) 'Xc=' num2str(round(sol.Xc*1e4)/10) ...
    'taup=' num2str(round(sol.tau_p*1e-4)/100) 'y=' num2str(abs(sol.y)*1e3)];
clear global mu ap bs k 

function F=calcFsym(tau,ILowerBoundery,z)
%Integral is done in non dimensional coordinats
%eq.6.3.49
global mu g ap

%tau_p and rp (Xc) are defined at begining
tau_p=1;
rp=1;%(Xc)
I=@(w) tau(w)./(w-z)./((-w).^(1-g));

%------Decomposition to sym and Asym
v=@(w) 1/2*( I(w)+ I(2*z-w));
%h=@(w) 1/2*( I(w)- I(2*z-w));
l=min(z-ILowerBoundery,0-z);
%I1=integral(v,z-l/2,z+l/2)+integral(I,-1,z-l/2)+integral(I,z+l/2,0,'RelTol',1e-3,'AbsTol',1e-5);
%I1=integral(I,-1,z-l/2)+2*integral(v,z-l/2,z)+my_integralSingular(I,1-g,z+l/2,0);
I1=integral(I,ILowerBoundery,z-l/2)+2*integral(v,z-l/2,z)+my_integralSingular(I,1-g,z+l/2,0);

F=tau_p*rp^(g-1)*( I1-1i*pi*tau(z)/(-z)^(1-g) );%Second part in (...) needed for Principal value.
z=z-1E-22i;      %allways work at lower half plane!!!!!!!!
F=1i*sin(pi*g)*z.^(1-g)/(2*pi*mu*ap).*F;

function F=calcFNumerical(tau,ILowerBoundery,z)
%Integral is done in non dimentional coordinats
%eq.6.3.49
global mu g ap

%-----Numerical integration
tau_p=1;
rp=1;%(Xc)

I=@(w) tau(w)./(w-z)./((-w).^(1-g));
%F=tau_p*rp^(g-1)*integral(I,ILowerBoundery,0,'RelTol',1e-2,'AbsTol',1e-25);%---without 'Tol' doesn't work for g<0.34 .
F=tau_p*rp^(g-1)*my_integralSingular(I,1-g,ILowerBoundery,0);
%F=tau_p*integral(I,-1,0)*rp^(g-1);
z=z-1E-22i;      %allways work at lower half plane!!!!!!!!
F=1i*sin(pi*g)*z.^(1-g)/(2*pi*mu*ap).*F;

function [Sxx,Syy,Sxy]=calcStress(F,f)
%eq6.3.31-6.3.33
global mu ap bs k
Sxy=mu*(1i*ap*(F-conj(F))+(bs^2-1)*f);
Syy=mu/2*(  (bs^2-1)*( F+conj(F) )+4*bs*f  );
Sxx=mu/2/k^2*(  (1-(1-2*k^2)*ap^2)*(F+conj(F))-4*k^2*bs*f );

%  Syy=mu/2*(  (bs^2-1)*( F+conj(F) ));%+4*bs*f  );
%  Sxx=mu/2/k^2*(  (1-(1-2*k^2)*ap^2)*(F+conj(F)));%-4*k^2*bs*f );
%
%  Syy=mu/2*( 4*bs*f  );
%  Sxx=mu/2/k^2*(-4*k^2*bs*f );


