function sol=CrackSolutionStokes(v,order,h)
%According to the Efim Brener's solution with stokes friction at the crack tail.
%The strains are calculated from the displacement field.
%[v]=[m/s]
%[h]=m

%--------
dx=1E-5;
sol.x=-70E-3:dx:70E-3;

dx=dx/10;  
x=-70E-3:dx:70E-3;%[m] distance from the crack tip ,theta=0;

[r,theta]=calc_r_theta(x,h);    
[ux,uy,n]=calc_u(v,r,theta,order); 
Uxx=diff(ux,1,2)/dx;
Uxy1=diff(uy,1,2)/dx;
clear ux uy

x=x+dx/2;
x=x(1:end-1);
[r,theta]=calc_r_theta(x,h-dx/2);    
[ux(1,:),uy(1,:)]=calc_u(v,r,theta,order); 
[r,theta]=calc_r_theta(x,h+dx/2);    
[ux(2,:),uy(2,:)]=calc_u(v,r,theta,order); 
Uyy=diff(uy,1,1)/dx;
Uxy2=diff(ux,1,1)/dx;

Uxy=1/2*(Uxy1+Uxy2);

%-----interpoltate the solution to coordinates as used in CrackSolutionForh.m  
sol.Uxy=interp1(x,Uxy,sol.x);
sol.Uxx=interp1(x,Uxx,sol.x);
sol.Uyy=interp1(x,Uyy,sol.x);
sol.n=n;
sol.v=v;




function [r,theta]=calc_r_theta(x,h)
%------find theta
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
%---------find r
r=(h^2+x.^2).^0.5;

function [ux,uy,n]=calc_u(v,r,theta,order) 

Cd=2700;
Cs=1345;
nu=(2*(Cs/Cd).^2-1)/2/((Cs/Cd).^2-1);%plane stress
ro=1.17E3;%[Kg/m^3];
E=Cs^2*2*ro*(1+nu);%[Pa]
mu=0.5*E/(1+nu);
G=1;%[J/m^2]

alpha_d=(1-(v/Cd).^2).^0.5;
alpha_s=(1-(v/Cs).^2).^0.5;
gamma_d=(1-(v*sin(theta)/Cd).^2).^0.5;
gamma_s=(1-(v*sin(theta)/Cs).^2).^0.5;
r_d=r.*gamma_d;
r_s=r.*gamma_s;

if (theta>0)
    theta_d=atan(alpha_d.*tan(theta));
    index=find(theta_d<0);%matlab gives -pi/2<theta<pi/2, I want 0<theta<pi
    theta_d(index)=theta_d(index)+pi;
    
    theta_s=atan(alpha_s.*tan(theta));
    index=find(theta_s<0);%matlab gives -pi/2<theta<pi/2, I want 0<theta<pi
    theta_s(index)=theta_s(index)+pi;
else
    theta_d=atan(alpha_d.*tan(theta));
    index=find(theta_d>0);%matlab gives -pi/2<theta<pi/2, I want -pi<theta<0
    theta_d(index)=theta_d(index)-pi;
    
    theta_s=atan(alpha_s.*tan(theta));
    index=find(theta_s>0);%matlab gives -pi/2<theta<pi/2, I want -pi<theta<0
    theta_s(index)=theta_s(index)-pi;
    
end

D=4*alpha_d.*alpha_s-(1+alpha_s.^2).^2;
A=alpha_s.*v.^2./D/(1-nu)/Cs^2;
K=(G*E/(1-nu^2)./A).^0.5; %plane strain

%-------calc n (\lambda)
%alpha=1;
%S=mu/alpha;
%S=Cs*10;
%n=1/pi*atan(S*D/2./v./alpha_s./(alpha_s.^2-1))+1;
l=0.78;
n=l+order;
%------calc displacement field
ux=2*K/(mu.*D*(2*pi)^0.5).*(2*alpha_s.*r_d.^n.*sin(n*theta_d)-alpha_s.*(1+alpha_s.^2)*r_s.^n.*sin(n*theta_s));       
uy=2*K/(mu.*D*(2*pi)^0.5).*(2*alpha_s.*alpha_d.*r_d.^n.*cos(n*theta_d)-(1+alpha_s.^2)*r_s.^n.*cos(n*theta_s));
