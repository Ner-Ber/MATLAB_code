function [tau_0 D]=Crack_SupershearSelfSimilarCalcTau0(a,b,k,mu,Gamma,d)

%Based on Broberg 1994.
ap=@(b)(1-b.^2).^0.5;
bs=@(b)((b/k).^2-1).^0.5;
g=@(b) 1/pi*atan(4*ap(b).*bs(b)./(bs(b).^2-1).^2);%g stands for gamma. Note that eq. 34 (book eq6.9.111) is for 1/b

D=@(w)(1+w);%Linear Cohesive zone model

%---------Calc wd (eq.67 book-eq.6.3.70)
%Linear cohesive zone 6.3.70 has two internal integrals,named here by I1  and I2
% for general cohesive zone this hould be an integral
I1=1/g(b);

%calculate I2
dz=1E-3;
z=0:dz:1-dz; % the integral of I is singular at z=1 but appears to cancel with the prefactord(-z) -> finite limit;

I2=zeros(1,length(z));

for l=1:length(z)
    I2_integrand=@(w) 1./(w-z(l))./((w).^(1-g(b)));
    I2(l)=integral(I2_integrand,1,1E7);
end

wd_integrand=D(-z).*z.^(1-g(b)).*(I1+D(-z).*I2);
wd=dz*trapz(wd_integrand);

%---------Calc N eq.47 (A.9)

%--------Calc Io and M
I=@(w) 2*(w.*g(1./w)-1/b*g(b))./(w.^2-1/b^2);
Io=integral(I,1,1/b,'RelTol',1e-6,'AbsTol',1e-3);%eq.41 (book eq.6.9.118)
M=sin(pi*g(b)).*exp(-Io).*b/2.*((1./b-1)./(1./b+1)).^g(b); %eq.46
%
ds=1e-5;
s_Vec=1+ds:ds:1/b;

for l=1:length(s_Vec)
    s=s_Vec(l);
    %-Io is a principle value integral.To solve it I decomposed to
    %symetric and asymetric parts. Asymetric part cancels out.Symetric
    %part is no longer singular.
    I=@(w) 2*(w.*g(1./w)-1/b*g(b))./(w.^2-s^2);
    v=@(w) 1/2*( I(w)+ I(2*s-w));
    %h=@(w) 1/2*( I(w)- I(2*s-w));
    
    if (s<=0.5*(s_Vec(1)+1/b))
        Io=integral(v,1,2*s-1);%eq.6.9.118
        Io=Io+integral(I,2*s-1,1/b);
    else
        Io=integral(v,2*s-1/b,1/b);%eq.6.9.118
        Io=Io+integral(I,1,2*s-1/b);
    end
    
    M_s(l)=sin(pi*g(1/s))*exp(-Io)/(1/b+s)*((s-1)/(s+1))^(g(b)/b/s); %eq.46
end

s=s_Vec;
N_integrand=( M*(2/b./(1/b-s)).^g(b)-M_s.*((1/b+s)./(1/b-s)).^(g(b)/b./s) ) ./(1/b-s);
N=ds*trapz(N_integrand)+M/g(b)*(2/b/(1/b-1))^g(b);%eq.A9
clear M_s;

f1=b.^2.*sin(pi*g(b))./(4*k^2*(1-b.^2).^0.5); %eq.54 very similar to YII
G1=f1.*sin(pi*g(b))*2/pi.*wd;
A=1./g(b)./(g(b)+1);%the integral in eq.70 is solved easily for linear cohesive zone
G2=pi*A.^-1.* b.*exp(-Io)./( 2.^(1-g(b)).*g(b).*N).*((1-b)./(1+b)).^g(b);
D=-pi/N;%eq.47 In the book 6.9.122 tau_0 doesn't apear. it comes with H.

%G0=pi*tau_0^2*a/mu/4/(1-k^2);%eq.74
%G=tau_0^2*a/mu*(d/a).^(1-2*g(b)).*(G1.*(G2).^2); %eq.68
tau_0=(1/Gamma*a/mu*(d/a).^(1-2*g(b)).*(G1.*(G2).^2))^-0.5; %eq.68


