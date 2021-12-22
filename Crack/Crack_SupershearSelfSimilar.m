function sol=Crack_SupershearSelfSimilar(v,a,y)
% b[m/s]
%Based on Broberg 1994.

if (y<0)
        y=-y; %to make sure that y >0
end

%--------dimensions
clear global b mu g ap bs k E nu PlaneStrain
global E nu b g PlaneStrain

PlaneStrain=true;
% ro=1.17E3;%[Kg/m^3];
% Cd=2700;
% Cs=1345;

%David simulation
ro=1130;
Cd=2017;
Cs=916;

mu=Cs^2*ro;
nu=(2*(Cs/Cd)^2-1)/2/((Cs/Cd)^2-1);%plane strain
E=Cs^2*2*ro*(1+nu);%[Pa]

if(PlaneStrain~=true)
    %-------Plain Stress
    Cd=(E/ro/(1-nu^2))^0.5; %approx. =2340m/s;
end

k=(Cs/Cd);
b=v/Cd;

%-----Define Gamma and Xc at2^0.5Cs
%Gamma=1; %[J/m^2]
%Xc0=5E-3;

%Davids Simulation
Gamma=3.4427;
Xc0=3.8e-3;

Xc_n0=Crack_SupershearCalcXc(2^0.5*Cs/Cd,Cs/Cd);
tau_d=(Xc_n0/Xc0*Gamma*mu)^0.5;
 Xc_n=Crack_SupershearCalcXc(v/Cd,Cs/Cd);
 Xc=Gamma*mu/tau_d^2*Xc_n;%Should be in [m]
[tau_0 D]=Crack_SupershearSelfSimilarCalcTau0(a,b,k,mu,Gamma,Xc);
%---------

ap=@(b) (1-b.^2).^0.5;
bs=@(b) ((b/k).^2-1).^0.5;
g=@(b) 1/pi*atan(4*ap(b).*bs(b)./(bs(b).^2-1).^2);%g stands for gamma. Note that eq. 34 (book eq6.9.111) is for 1/b
R=@(z) 4*k^3*(1-z.^2).^0.5.*(k^2-z.^2).^0.5-(z.^2-2*k^2).^2;

%I_int=@(w,z) 2*w.*g(1./w)./(w.^2-z^2);%eq.6.9.115
%I=@(z) integral(@(w) I_int(w,z),1,1/b);

H0=@(z) exp(-I(z));%eq.6.9.114
H=@(z) D*H0(z)./(b^-2-z.^2);%eq.6.9.120
A=@(z) tau_0*2*k^4*(k^-2 - z.^2).^0.5./(mu*R(1/z))*H(z);%Factor of q^-2 is missing.
C=@(z)-tau_0*k^4*(k^-2 - 2*z.^2)./(mu*z*R(1/z))*H(z); %Factor of q^-2 is missing.

%----------Following the book y>0 (y is allready has been defined)
r=@(x,y) (x.^2+y.^2).^0.5;
z_p=@(x,s,k) r(x,y).^-2 .*(x*s+1i*y/k*(k^2*s^2-r(x,y).^2).^0.5);%book 6.9.10
g_p=@(x,s,k) r(x,y).^-2 .*(x+1i*y*k*s*(k^2*s^2-r(x,y).^2 ).^-0.5);%book 6.9.12

t=a/b;
%t=xf;%tau
%x_vec=0.1*t*k:1e-4:t;
%x_vec=t-200e-3:1e-4:t+60e-3;
%x_vec=0.1*t*k:1e-5:0.9*t*b;
x_vec=t*k:0.5e-5:0.99*t;
%--calculate the derivitves of the potentials (S and F)

%xi_m=@(x,y,t) 1./r(x,y).*(  -x*t./r(x,y)-( (x*t./r(x,y)).^2-t^2+(y/k).^2 ).^0.5 );
xi_p=@(x,y,t) 1./r(x,y).*(  -x*t./r(x,y)+( (x*t./r(x,y)).^2-t^2+(y/k).^2 ).^0.5 );


%for check
%f=@(xi,x) xi*x-(k^-2-xi.^2).^0.5*y;
%l=1;
% for j=1:length(x_vec)
%     x=x_vec(j);
%     if (x+(k^-2-1)^0.5*y<t)&&(t<x^2/k/r(x,y)+(k^-2-(x/k/r(x,y))^2)^0.5*y)
%        xx(l)=x;
%        l=l+1;
%     end
%
%     ee=-x/k/r(x,y):1e-5:-1;
%     [~,index]=min(abs(f(ee,x)+t));
%     xi_m2=
% end

for j=1:length(x_vec)
    x=x_vec(j);
    
    if(t-r(x,y)/k>0) %distances where Cs have already arrived.
        S_xxx(j)=-1/pi*imag(g_p(x,-t,k)*C(z_p(x,-t,k)));%book 6.9.15
        
        ss=-z_p(x,-t,k)^-1*(k^-2-z_p(x,-t,k)^2)^0.5;
        S_xyx(j)=-1/pi*imag(g_p(x,-t,k)*ss*C(z_p(x,-t,k)));%book 6.9.15
        
        ss=z_p(x,-t,k)^-2*(k^-2-z_p(x,-t,k)^2);%pre factor due to derivitive
        S_yyx(j)=-1/pi*imag(g_p(x,-t,k)*ss*C(z_p(x,-t,k)));%book 6.9.15
        
    else
        S_xxx(j)=0;
        S_yyx(j)=0;
        S_xyx(j)=0;
    end
    
    if (x+(k^-2-1)^0.5*y<t)&&(t<x^2/k/r(x,y)+(k^-2-(x/k/r(x,y))^2)^0.5*y) 
        
        S_xxx(j)=S_xxx(j)--1/pi*imag(C(xi_p(x,y,t)));%book 6.9.15
        
        ss=-xi_p(x,y,t)^-1*(k^-2-xi_p(x,y,t)^2)^0.5;
        S_xyx(j)=S_xyx(j)--1/pi*imag(ss*C(xi_p(x,y,t)));%book 6.9.15
        
        ss=xi_p(x,y,t)^-2*(k^-2-xi_p(x,y,t)^2);%pre factor due to derivitive
        S_yyx(j)=S_yyx(j)--1/pi*imag(ss*C(xi_p(x,y,t)));%book 6.9.15
        
    end
   
    
    if(t-r(x,y)>0) %distances where Cp have already arrived.
        
        F_xxx(j)=-1/pi*imag(g_p(x,-t,1)*A(z_p(x,-t,1)));%book 6.9.15
        
        ss=-z_p(x,-t,1)^-1*(1-z_p(x,-t,1)^2)^0.5;
        F_xyx(j)=-1/pi*imag(g_p(x,-t,1)*ss*A(z_p(x,-t,1)));%book 6.9.15
        
        ss=z_p(x,-t,1)^-2*(1-z_p(x,-t,1)^2);%pre factor due to derivitive
        F_yyx(j)=-1/pi*imag(g_p(x,-t,1)*ss*A(z_p(x,-t,1)));%book 6.9.15
        
    else
        F_xxx(j)=0;
        F_xyx(j)=0;
        F_yyx(j)=0;
        
    end
    
end

% for j=1:length(x_vec) txy(j)=-1/pi/x_vec(j)*imag(H(-t/x_vec(j))) ; end
% sol.txy=txy;

Sxy_x=mu*(2*F_xyx+S_yyx-S_xxx);
Syy_x=mu/k^2*( (1-2*k^2)*F_xxx+F_yyx-2*k^2*S_xyx );
Sxx_x=mu/k^2*(F_xxx+(1-2*k^2)*F_yyx+2*k^2*S_xyx );

%-------Works for both boundary conditions with appropriate k 
 sol.Uxy_x=1/(2*mu)*Sxy_x;
 sol.Uxx_x=1/(2*mu)/(2*(1-k^2))*(Sxx_x-(1-2*k^2)*Syy_x);
 sol.Uyy_x=1/(2*mu)/(2*(1-k^2))*(Syy_x-(1-2*k^2)*Sxx_x);

 %---------Integrate (doesn't converg!)
dx=x_vec(2)-x_vec(1);
Sxy=cumsum(Sxy_x)*(x_vec(2)-x_vec(1))*dx;
Sxy=Sxy(end)-Sxy;

Syy=cumsum(Syy_x)*(x_vec(2)-x_vec(1))*dx;
Syy=Syy(end)-Syy;

Sxx=cumsum(Sxx_x)*(x_vec(2)-x_vec(1))*dx;
Sxx=Sxx(end)-Sxx;

%-------Works for both boundary conditions with appropriate k 
 sol.Uxy=1/(2*mu)*Sxy;
 sol.Uxx=1/(2*mu)/(2*(1-k^2))*(Sxx-(1-2*k^2)*Syy);
 sol.Uyy=1/(2*mu)/(2*(1-k^2))*(Syy-(1-2*k^2)*Sxx);

 

sol.v=v;
sol.Cd=Cd;
sol.Cs=Cs;
sol.mu=mu;
sol.x=x_vec;
sol.a=a;
sol.y=y;
sol.Sxy_x=Sxy_x;
sol.Syy_x=Syy_x;
sol.Sxx_x=Sxx_x;
sol.Sxy=Sxy;
sol.Sxx=Sxx;
sol.Syy=Syy;
sol.Xc=Xc;
clear global b g;

function I=I(z)
global b g;

if (imag(z)~=0 )
    
    I_int=@(w) 2*w.*g(1./w)./(w.^2-z^2);%eq.6.9.115
    I=integral( I_int ,1,1/b);
else
 
    if ~((1<abs(z))&&(abs(z)<1/b))
        I_int=@(w) 2*w.*g(1./w)./(w.^2-z^2);%eq.6.9.115
        I=integral(I_int,1,1/b);
        
    else
       
        %-Io is a principle value integral.To solve it I decomposed to
        %symetric and asymetric parts. Asymetric part cancels out.Symetric
        %part is no longer singular.
        z=abs(z);
        
        I_int=@(w) 2*w.*g(1./w)./(w.^2-z^2);
        v=@(w) 1/2*( I_int(w)+ I_int(2*z-w));
        %h=@(w) 1/2*( I(w)- I(2*s-w));
        
        if (z<=0.5*(1+1/b))
            I=integral(v,1,2*z-1);%eq.6.9.118
            I=I+integral(I_int,2*z-1,1/b);
        else
            I=integral(v,2*z-1/b,1/b);%eq.6.9.118
            I=I+integral(I_int,1,2*z-1/b);
        end

        I=I+pi*1i*g(1./z) ;
        
    end
end



