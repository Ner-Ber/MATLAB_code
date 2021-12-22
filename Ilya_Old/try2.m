g=sym(1/5);
z=-0.4;

syms w;
tau(w)=(1+w); %linear weakening cohesive zone in Normlized Units(w=z/Xc)
I=tau/(w-z)/((-w)^(1-g));
Io=int(I,w,-1,0,'PrincipalValue', true);
Io=eval(Io)
clear I tau

g=eval(g);
tau=@(w) 1+w;
I=@(w) tau(w)./(w-z)./((-w).^(1-g));
v=@(w) 1/2*( I(w)+ I(2*z-w));
%h=@(w) 1/2*( I(w)- I(2*z-w));
l=min(z--1,0-z);
Io=integral(v,z-l/2,z+l/2)+integral(I,-1,z-l/2)+integral(I,z+l/2,0,'RelTol',1e-3,'AbsTol',1e-5)
Io=integral(v,z-l/2,z+l/2)+integral(I,-1,z-l/2)+quad(I,z+l/2,0)
Io=integral(v,z-l/2,z+l/2)+integral(I,-1,z-l/2)+my_integralSingular(I,1-g,z+l/2,0)
clear