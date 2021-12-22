function sol=Crack_Supershear(v,y)

% !!!!!!!!Lower half plane-Use only y<0!!!!!!!!
%Broberg p348

%----------Calculations done only on the lower half plane,y<0.
%If y>0 , calculations done for y<0 and then Sxx and Syy multiplied by -1.

if (y>0)
    LowerHalfPlane=false;
    y=-y; %to make sure that y <0
else
    LowerHalfPlane=true;
end


clear global mu g ap bs k E nu rp PlaneStrain;
global mu g ap bs k E nu rp PlaneStrain
%warning('off');
PlaneStrain=false;

Cd=2700;
Cs=1345;
ro=1.17E3;%[Kg/m^2];

nu=(2*(Cs/Cd).^2-1)/2/((Cs/Cd).^2-1);%
E=Cs^2*2*ro*(1+nu);%[Pa]

if(PlaneStrain~=true)
    %-------Plain Stress
    Cd=(E/ro/(1-nu^2))^0.5; %approx. =2340m/s;
end

%-----------
mu=Cs^2*ro;
%k=((1-2*nu)/2/(1-nu))^0.5;%plane strain
%k=((1-nu)/2)^0.5;%plane stress
k=(Cs/Cd);
b=v/Cd;
%b=v/Cd;%p332
ap=(1-b^2)^0.5;
bs=((b/k)^2-1)^0.5;
alpha=bs/ap;
g=1/pi*atan(4*ap*bs/(bs^2-1)^2);

%-----Define Gamma and Xc at2^0.5Cs
Gamma=1; %[J/m^2]
Xc0=5E-3;
Xc_n0=Crack_SupershearCalcXc(2^0.5*Cs/Cd,Cs/Cd);
tau_s=(Xc_n0/Xc0*Gamma*mu)^0.5;
Xc_n=Crack_SupershearCalcXc(v/Cd,Cs/Cd);
rp=Gamma*mu/tau_s^2*Xc_n;%Should be in [m]

%rp=1.3E-3;
%tau_s=1.05E6;
%---------------

dx=1E-6;
x=70E-3:-dx:-70E-3;%[m] distance from the crack tip ,theta=0;
y=ap*y; %normlized units eq.6.3.22

%------------Far field calculation
Ftmp=calcFarFieldF(x-alpha*y);
f=-(bs^2-1)/(4*bs)*(Ftmp+conj(Ftmp));
F=calcFarFieldF(x+1i*y);

f=f*tau_s;
F=F*tau_s;

[sol.Sxx,sol.Syy,sol.Sxy]=calcStress(F,f);

if(LowerHalfPlane==false)
    sol.Sxx=-sol.Sxx;
    sol.Syy=-sol.Syy;
    sol.y=y/ap*rp;
    y=-y;
end

%-------Works for both boundary conditions with appropriate k 
 sol.Uxy=1/(2*mu)*sol.Sxy;
 sol.Uxx=1/(2*mu)/(2*(1-k^2))*(sol.Sxx-(1-2*k^2)*sol.Syy);
 sol.Uyy=1/(2*mu)/(2*(1-k^2))*(sol.Syy-(1-2*k^2)*sol.Sxx);


sol.phi_f=asin(k/b)+pi/2;
sol.x=x;
sol.xMach=y*alpha;
sol.y=y/ap;
sol.v=v;
sol.g=g;
sol.Cs=Cs;
sol.Cd=Cd;
sol.Cr=NaN;

clear global mu g ap bs k E nu rp;

function F=calcFarFieldF(z)
%--Far field solutions. 6.3.55
global mu g ap rp
tau_s=1;
tau=@(w)(1+w);
I=@(w) tau(-w)./(w).^(1-g);
%A=sin(pi*g)/pi*tau_s*integral(I,0,1)*rp^g;

A=sin(pi*g)/pi*tau_s*integral(I,0,1,'RelTol',1e-2,'AbsTol',1e-25)*rp^g;
%A=5.2E4;
z=z-1E-22i;     %allways work at lower half plane!!!!!!!!
F=-1i*A./(2*mu*ap*z.^g);


function [Sxx,Syy,Sxy]=calcStress(F,f)
global mu ap bs k
Sxy=mu*(1i*ap*(F-conj(F))+(bs^2-1)*f);
Syy=mu/2*(  (bs^2-1)*(F+conj(F))+4*bs*f  );
Sxx=mu/2/k^2*(  (1-(1-2*k^2)*ap^2)*(F+conj(F))-4*k^2*bs*f );

