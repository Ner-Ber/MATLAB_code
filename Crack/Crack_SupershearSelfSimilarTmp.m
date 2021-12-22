function sol=Crack_SupershearSelfSimilarTmp(b)
% b[v/Cd]
%Based on Broberg 1994.

%--------dimensions
clear global mu g ap bs k E nu PlaneStrain
global E nu PlaneStrain

PlaneStrain=false;
ro=1.17E3;%[Kg/m^3];
Cd=2700;
Cs=1345;

mu=Cs^2*ro;
nu=(2*(Cs/Cd)^2-1)/2/((Cs/Cd)^2-1);%plane strain
E=Cs^2*2*ro*(1+nu);%[Pa]

if(PlaneStrain~=true)
    %-------Plain Stress
    Cd=(E/ro/(1-nu^2))^0.5; %approx. =2340m/s;
end

k=(Cs/Cd);
D=1;

ap=@(b)(1-b.^2).^0.5;
bs=@(b)((b/k).^2-1).^0.5;
g=@(b) 1/pi*atan(4*ap(b).*bs(b)./(bs(b).^2-1).^2);%g stands for gamma. Note that eq. 34 (book eq6.9.111) is for 1/b
%R=@(z) 4*k^3*(1-z.^2).^0.5.*(k^2-z.^2).^0.5-(z.^2-2*k^2).^2;

R=@(z) 4*z.^2.*(1-z.^2).^0.5.*(k^-2-z.^2).^0.5+(k^-2-2*z.^2).^2;%It is not exactlly as in the book

I_int=@(w,z) 2*w.*g(1./w)./(w.^2-z^2);
I=@(z) integral(@(w) I_int(w,z),1,1/b,'RelTol',1e-6,'AbsTol',1e-5);
H0=@(z) exp(-I(z));
H=@(z) D*H0(z)./(b^-2-z.^2);
%A=@(z)-2*k^4*(k^-2 - z.^2).^0.5 ./(mu*R(1/z))*H(z);%Factor of q^-2 is missing
%C=@(z)-k^4*(k^-2 - 2*z.^2)./(mu*z*R(1/z))*H(z); %Factor of q^-2 is missing

A=@(z) -2*z.^4*(k^-2 - z.^2).^0.5./mu/R(z)*H(z);%Factor of q^-2 is missing-my
C=@(z) z.^3*(k^-2 -2*z.^2)./mu/R(z)*H(z);%Factor of q^-2 is missing-my

%----------Following the book y>0
y=3.4e-3;
t=100e-3;

z_p=@(x,s,k) (x.^2+y.^2).^-1 .*(x*s+1i*y/k*(k^2*s^2-(x.^2+y.^2)).^0.5);%book 6.9.10
g_p=@(x,s,k) (x.^2+y.^2).^-1 .*(x+1i*y*k*s*(k^2*s^2-(x.^2+y.^2)).^-0.5);%book 6.9.12

x_vec=0.6*t*k:0.01e-3:t;%tau

for j=1:length(x_vec)
    x=x_vec(j);
    r=(x.^2+y.^2)^0.5;
    
    if(t-r/k>0)
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
    
    if(t-r>0)
        %I don't understand but with additional "-" it makes more sense 
        F_xxx(j)=-1/pi*imag(g_p(x,-t,1)*A(z_p(x,-t,1)));%book 6.9.15
        
        ss=-z_p(x,-t,1)^-1*(1-z_p(x,-t,1)^2)^0.5;%pre factor due to derivitive in y direction
        F_xyx(j)=-1/pi*imag(g_p(x,-t,1)*ss*A(z_p(x,-t,1)));%book 6.9.15
        
        ss=z_p(x,-t,1)^-2*(1-z_p(x,-t,1)^2);%pre factor due to derivitive in y direction
        F_yyx(j)=-1/pi*imag(g_p(x,-t,1)*ss*A(z_p(x,-t,1)));%book 6.9.15
    else
        F_xxx(j)=0;
        F_xyx(j)=0;
        F_yyx(j)=0;
        
    end
    
end


Sxy_x=mu*(2*F_xyx+S_yyx-S_xxx);
Syy_x=mu/k^2*( (1-2*k^2)*F_xxx+F_yyx-2*k^2*S_xyx );
Sxx_x=mu/k^2*(F_xxx+(1-2*k^2)*F_yyx+2*k^2*S_xyx );

Sxy=cumsum(Sxy_x)*(x_vec(2)-x_vec(1));
Sxy=Sxy(end)-Sxy;

Syy=cumsum(Syy_x)*(x_vec(2)-x_vec(1));
Syy=Syy(end)-Syy;

Sxx=cumsum(Sxx_x)*(x_vec(2)-x_vec(1));
Sxx=Sxx(end)-Sxx;

[sol.Uxx,sol.Uyy,sol.Uxy]=calcStrain(Sxx,Syy,Sxy);

sol.x=x_vec;
sol.Sxy_x=Sxy_x;
sol.Syy_x=Syy_x;
sol.Sxy=Sxy;
sol.Sxx=Sxx;
sol.Syy=Syy;

clear global E nu PlaneStrain

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

