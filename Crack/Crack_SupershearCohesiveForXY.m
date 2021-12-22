function sol=Crack_SupershearCohesiveForXY(g_in,x,y)
%The function calculates elastic fields of super shear rupture according to
%Broberg p348. Assumes linear cohesive zone  weakening. 
%Inseatd of specifying rupture velocity enter g_in (eq.6.3.42). g_in - should be symbolic. sym(1/3),sym(5/11)... 
%[y]=m,[x]=m can be vectores

%----------Calculations done only on the lower half plane,y<0.
%If y>0 , calculations done for y<0 and then Sxx and Syy multiplied by -1.

if (y>0)
    LowerHalfPlane=false;
    y=-y; %to make sure that y <0
else
    LowerHalfPlane=true;
end
%----------------
clear global mu g_sym ap bs k E nu PlaneStrain
global mu g_sym ap bs k E nu PlaneStrain
g_sym=g_in;
PlaneStrain=false;

%--------Elastic properties
ro=1.17E3;%[Kg/m^3];
Cd=2700;
Cs=1345;

%David simulation
%ro=1130;
%Cd=2017;
%Cs=916;

nu=(2*(Cs/Cd)^2-1)/2/((Cs/Cd)^2-1);%plane strain
E=Cs^2*2*ro*(1+nu);%[Pa]

if(PlaneStrain~=true)
    %-------Plain Stress
    Cd=(E/ro/(1-nu^2))^0.5; %approx. =2340m/s;
end

%--------calc v and g
%v=v*Cs;
% g=1/pi*atan(4*ap*bs/(bs^2-1)^2);
f=@(Cf) 4./( (Cf/Cs).^2-2 ).^2 .*( ( 1-(Cf/Cd).^2 ).*(  (Cf/Cs).^2-1) ).^0.5;
%g=1/pi*atan(f(v));
%g_sym=sym(1/2);
g=eval(g_sym);
%--Cf>2^0.5Cs
v_tmp=2^0.5*Cs:1E-4*Cs:Cd;
[~,index]=min(abs(tan(pi*g)-f(v_tmp)));
v=v_tmp(index);
%---------------------

k=(Cs/Cd);%general
%k=((1-2*nu)/2/(1-nu))^0.5;%plane strain
%k=((1-nu)/2)^0.5;%plane stress
mu=Cs^2*ro;
b=v/Cd;%p332
ap=(1-b.^2).^0.5;
bs=((b/k).^2-1).^0.5;
alpha=bs/ap;
phi_f=asin(k/b);
%------Frictional properties
%Unlike Sub-Rayleigh two frictional properties should be specified.

syms w;
 tau(w)=(1+w); %linear weakening cohesive zone in Normlized Units(w=z/Xc)
 ILowerBoundery=-1;%Lower boundary of the integration. cohesive zone size

%Symbolic integration doesn't work 
 % tau(w)=exp(w);
% ILowerBoundery=-1;

%------Define Gamma (fracture energy) and Xc
Gamma=1.1; %[J/m^2]
Xc=1.1E-3;%Xc=3E-3;
% 
% %Gamma=3.4427;
% %Xc=2.75E-3;

Xc_n=Crack_SupershearCalcXc(v/Cd,Cs/Cd);%eq6.3.69
tau_s=(Xc_n/Xc*Gamma*mu)^0.5;%peak shear stress
Xc_n0=Crack_SupershearCalcXc(2^0.5*Cs/Cd,Cs/Cd);
Xc0=Xc_n0*Gamma*mu/tau_s^2;

%------Define Gamma and tau_s
%  Gamma=1; %[J/m^2]
%  tau_s=1.05E6;%tau_s=mu/1000;%maximal shear stress
% 
% %David's simulation
% %Gamma=3.4427;
% %tau_s=1.5792e6;
% 
% Xc_n=Crack_SupershearCalcXc(v/Cd,Cs/Cd);
% Xc=Gamma*mu/tau_s^2*Xc_n;%Should be in [m]

%-----Define Gamma and Xc at2^0.5Cs
% Gamma=1.1; %[J/m^2]
% Xc0=3E-3;
% Xc_n0=Crack_SupershearCalcXc(2^0.5*Cs/Cd,Cs/Cd);
% tau_s=(Xc_n0/Xc0*Gamma*mu)^0.5;
% Xc_n=Crack_SupershearCalcXc(v/Cd,Cs/Cd);
% Xc=Gamma*mu/tau_s^2*Xc_n;%Should be in [m]

%----Define x and y vectors - comment this part if you use x,y as input
xStart=-25e-3;
l=-xStart*sin(phi_f);
xEnd=xStart+l*sin(phi_f)+5e-3;
x=linspace(xStart,xEnd,600);
y=-(x-xStart)*tan(phi_f)^-1;

%----x and y should be normlized by Xc
y=ap*y/Xc; %normlized units eq.6.3.22
x=x/Xc;
%x=x+alpha*y;%to make symbolic calculation faster
 
%------------Cohesive zone calculation
F=zeros(1,length(x));
f=zeros(1,length(x));
for j=1:length(x)
    
    %--------calc F_tmp for later calculation of f
    if ((x(j)-alpha*y(j)<0 ))
        if (x(j)-alpha*y(j)>=ILowerBoundery)
            (x(j)-alpha*y(j))  %to follow the calculation
            Ftmp=calcFsym(tau,ILowerBoundery,x(j)-alpha*y(j));
            f(j)=-(bs^2-1)/(4*bs)*(Ftmp+conj(Ftmp));%eq.6.3.37
        else
            Ftmp=calcFNumerical(tau,ILowerBoundery,x(j)-alpha*y(j));
            f(j)=-(bs^2-1)/(4*bs)*(Ftmp+conj(Ftmp));
        end
    else
        f(j)=0;
        Ftmp=0;
    end
    
    %--------calc F eq.6.3.49 
    %Symbolic calculation is long. When it is not necessary numerical
    %integration is used.
    if (y(j)~=0)
        F(j)=calcFNumerical(tau,ILowerBoundery,x(j)+1i*y(j));
    elseif (x(j)>0)
        F(j)=calcFNumerical(tau,ILowerBoundery,x(j)+1i*y(j));
    else
        F(j)=Ftmp;
    end
    
end

f=f*tau_s;
F=F*tau_s;

[soltmp.Sxx,soltmp.Syy,soltmp.Sxy]=calcStress(F,f);
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
% 
% sol.x=xStart:dx*Xc:xEnd;
% sol.y=y/ap*Xc;
% sol.Sxx=spline(soltmp.x,soltmp.Sxx,sol.x);
% sol.Syy=spline(soltmp.x,soltmp.Syy,sol.x);
% sol.Sxy=spline(soltmp.x,soltmp.Sxy,sol.x);

sol=soltmp;
%----------------

[sol.Uxx,sol.Uyy,sol.Uxy]=calcStrain(sol.Sxx,sol.Syy,sol.Sxy);

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
sol.tau_s=tau_s;
sol.smtT=ceil(1E-6/ (abs( mean( diff(sol.x))) / sol.v) );

clear global mu g_sym ap bs k E nu PlaneStrain

function F=calcFsym(tau,ILowerBoundery,z)
%Integral is done in non dimensional coordinats
%eq.6.3.49
global mu g_sym ap

g=g_sym;
tau_s=1;
rp=1;%(Xc)
syms w;
% tau(w)=(1+w); %linear weakening cohesive zone
I=tau/(w-z)/((-w)^(1-g));
I=int(I,w,ILowerBoundery,0,'PrincipalValue', true);
I=eval(I);
g=eval(g);
F=tau_s*rp^(g-1)*( I-1i*pi*eval(tau(z))/(-z)^(1-g) );%Second part in (...) needed for Principal value.
z=z-1E-22i;      %allways work at lower half plane!!!!!!!!
F=1i*sin(pi*g)*z.^(1-g)/(2*pi*mu*ap).*F;

function F=calcFNumerical(tau,ILowerBoundery,z)
%Integral is done in non dimentional coordinats
%eq.6.3.49
global mu g_sym ap

g=eval(g_sym);
%-----Numerical integration
tau_s=1;
rp=1;%(Xc)
tau=matlabFunction(tau); 
%tau=@(w)(1+w);

I=@(w) tau(w)./(w-z)./((-w).^(1-g));
F=tau_s*rp^(g-1)*integral(I,ILowerBoundery,0,'RelTol',1e-2,'AbsTol',1e-25);%---without 'Tol' doesn't work for g<0.34 .
%F=tau_s*integral(I,-1,0)*rp^(g-1);
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

 
function [Uxx,Uyy,Uxy]=calcStrain(Sxx,Syy,Sxy)
global E nu PlaneStrain
Uxy=(1+nu)/E*Sxy;

if (PlaneStrain==true)
    %planes strain
    Uxx=1/E*(1+nu)*(Sxx*(1-nu)-nu*Syy);
    Uyy=1/E*(1+nu)*(Syy*(1-nu)-nu*Sxx);
else
    %planes stress
    Uxx=1/E*(Sxx-nu*Syy);
    Uyy=1/E*(Syy-nu*Sxx);
end
