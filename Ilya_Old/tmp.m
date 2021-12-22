function Ks=tmp(x,tau_0,a,g)

%Ks=tau_0*a.^g/g; %homogeneous case

tau_0=@(s) interp1(x,tau_0,s);
Ks=a*0;
for j=1:length(g)
    I=@(s) tau_0(a-s)./(s).^(1-g(j));
    %Ks(j)= my_integralSingular(I,1-g(j),0,a);
    Ks(j)= integral(I,0,a);
end