function sol=CrackSolution_SelfSimilar(v,xf,y)
% v [m/s]
%Based on Broberg book
% First the derivitive with respect to x of all the fields are calculated.
% i.e Sxx_x,Syy_x,Sxy_x. Then simple integration (cumsum) is done. One
% should be very carefull with that.

%--------dimensions
clear global mu g ap bs k E nu PlaneStrain
global E nu PlaneStrain

PlaneStrain=true;
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
D=-1;%otherwise all components are negative. Need to calculate
b=v/Cd;

R=@(z) 4*k^3*(1-z.^2).^0.5.*(k^2-z.^2).^0.5-(z.^2-2*k^2).^2;
H=@(z) D*z.^4.*R(1./z)./(k^-2-z.^2).^0.5./(b^-2-z.^2).^1.5;%eq.6.9.39+note on p.414
A=@(z) 2*k^4*(k^-2 - z.^2).^0.5./(R(1./z)).*H(z);%Factor of q^-2 is missing
C=@(z)-k^4*(k^-2 - 2*z.^2)./(z.*R(1./z)).*H(z); %Factor of q^-2 is missing


%--mode I
%H=@(z) D*z^4*R(1/z)/(1-z.^2).^0.5/(b^-2-z^2).^1.5;
% A=@(z)-k^4*(k^-2 - 2*z.^2)./(mu*z*R(1/z))*H(z);%6.9.34 Factor of q^-2 is missing
% C=@(z)-2*k^4*(1 -z.^2).^0.5./(mu*R(1/z))*H(z); %6.9.35 Factor of q^-2 is missing

%----------Following the book y>0
% y=3.4e-3;
% xf=30e-3;

z_p=@(x,s,k) (x.^2+y.^2).^-1 .*(x*s+1i*y/k*(k^2*s^2-(x.^2+y.^2)).^0.5);%book 6.9.10
g_p=@(x,s,k) (x.^2+y.^2).^-1 .*(x+1i*y*k*s*(k^2*s^2-(x.^2+y.^2)).^-0.5);%book 6.9.12

t=xf/v*Cd;
x_vec=0.1*t*k:1e-5:t;

% %another option is integrate by using integral and not cumsum
%  S_xxx1=@(x) heaviside(t-(x.^2+y.^2).^0.5/k).*-1/pi.*imag(g_p(x,-t,k).*C(z_p(x,-t,k)));%book 6.9.15
%  S_xx=@(x) integral(S_xxx1,x_vec(1),x);
 
 
for j=1:length(x_vec)
    x=x_vec(j);
    r=(x.^2+y.^2)^0.5;
    
    if(t-r/k>0)
        S_xxx(j)=-1/pi*imag(g_p(x,-t,k)*C(z_p(x,-t,k)));%book 6.9.15
%         S_xxV(j)=S_xx(x);
        ss=-z_p(x,-t,k)^-1*(k^-2-z_p(x,-t,k)^2)^0.5;%pre factor due to derivitive
        S_xyx(j)=-1/pi*imag(g_p(x,-t,k)*ss*C(z_p(x,-t,k)));%book 6.9.15
        
        ss=z_p(x,-t,k)^-2*(k^-2-z_p(x,-t,k)^2);%pre factor due to derivitive
        S_yyx(j)=-1/pi*imag(g_p(x,-t,k)*ss*C(z_p(x,-t,k)));%book 6.9.15
        
    else
        
%         S_xxV(j)=S_xxV(j-1);
        S_xxx(j)=0;
        S_yyx(j)=0;
        S_xyx(j)=0;
    end
    
    if(t-r>0)
        F_xxx(j)=-1/pi*imag(g_p(x,-t,1)*A(z_p(x,-t,1)));%book 6.9.15
        
        ss=-z_p(x,-t,1)^-1*(1-z_p(x,-t,1)^2)^0.5;%pre factor due to derivitive
        F_xyx(j)=-1/pi*imag(g_p(x,-t,1)*ss*A(z_p(x,-t,1)));%book 6.9.15
        
        ss=z_p(x,-t,1)^-2*(1-z_p(x,-t,1)^2);%pre factor due to derivitive
        F_yyx(j)=-1/pi*imag(g_p(x,-t,1)*ss*A(z_p(x,-t,1)));%book 6.9.15
    else
        F_xxx(j)=0;
        F_xyx(j)=0;
        F_yyx(j)=0;
        
    end
    
end

S_xxV2=cumsum(S_xxx)*(x_vec(2)-x_vec(1));

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
sol.xf=xf;
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

