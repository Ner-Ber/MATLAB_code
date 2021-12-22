function [N,I]=ContactImaging(theta,lambda)
%incident angle (deg)
%incident wave length [sigma]
%theta=80;

sigma=1;
Phi= @(z) 1/(2*pi)^0.5/sigma*exp(-0.5*(z/sigma).^2); %distribution of heights.
lambda=sigma*lambda;%10,1
h_vec=linspace(1e-1*sigma,5*sigma,100);

for j=1:length(h_vec)
    
    h=h_vec(j);
    %Integrand = @(z) (z-h).*Phi(z); %Total area from elatic response
    Integrand = @(z) Phi(z);
    N(j)=integral(Integrand,h,1e4);
    
    %accurate FTIR
    Integrand = @(z) Phi(z).*FTIR(theta,(h-z)/lambda);

    %Exponential approximation
%     n1=1.5;
%     n2=1;
%     l=lambda/(4*pi)./(n1^2*sind(theta).^2-n2^2).^0.5;
%     Integrand = @(z) Phi(z).*exp(-(h-z)/l);
    
    I(j)=integral(Integrand,-1e4,h);
   
    Integrand = @(z) Phi(z); 
    I(j)=I(j)+integral(Integrand,h,1e4);
end

figure(20);
plot(N,I,'.-');
my_legend_add(['theta=' num2str(theta) 'lambda=' num2str(lambda) ]);
