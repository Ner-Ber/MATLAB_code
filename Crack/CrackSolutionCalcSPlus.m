function SPlus=CrackSolutionCalcSPlus(h,a,b)
%[h]=1/(m/s)
%Following freund Eq. 6.4.18 -Calculating S_+(h) -> S_+(h)=1/S_-^0(h) p.349

alpha=@(z) ( a^2-z.^2 ).^0.5;
beta=@(z) (b^2-z.^2).^0.5;
V=@(z) 4*z.^2.*abs(beta(z)).*abs(alpha(z)) ./ (2*z.^2-b^2).^2;
aPlus=a;%p344
bPlus=b;
SPlus=h*0; %allocation

for j=1:length(h)
    f=@(z) -1/pi*atan(V(z))./(z-h(j)); %Integrand of eq.6.4.18
    %z=linspace(aPlus,bPlus,1000);
    S0Minus=exp(integral(f,aPlus,bPlus));
    SPlus(j)=1/S0Minus;
end

% %Following freund Eq. 6.4.18 -Calculating S_+(h) -> S_+(h)=1/S_-^0(h) p.349
%
% alpha=@(z) ( a^2-z.^2+a^2*z.^2/h^2-2*a^2*z/h ).^0.5;
% beta=@(z) (b^2-z.^2+b^2*z.^2/h^2-2*b^2*z/h).^0.5;
% V=@(z) 4*z.^2.*abs(beta(z)).*abs(alpha(z)) ./ (2*z.^2-b^2).^2;
%
% f=@(z) -1/pi*atan(V(z))./(z-h); %Integrand of eq.6.4.18
%
% aPlus=a/(1+a/h);%p344
% bPlus=b/(1+b/h);
% z=linspace(aPlus,bPlus,1000);
% S0Minus=exp(integral(f,aPlus,bPlus));
% SPlus=1/S0Minus;

% %---------------------S_+ Not working
% alpha=@(z) ( a^2-z.^2+a^2*z.^2/h^2-2*a^2*z/h ).^0.5;
% beta=@(z) (b^2-z.^2+b^2*z.^2/h^2-2*b^2*z/h).^0.5;
% V=@(z) 4*z.^2.*beta(z).*abs(alpha(z)) ./ (2*z.^2-b^2-b^2*z.^2/h^2-2*b^2*z/h).^2;
%
% f=@(z) atan(V(z))./(z+h); %Integrand of eq.6.4.18
%
% aMin=a/(1-a/h);%p344
% bMin=b/(1-b/h);
%
% SPlus=exp(-1/pi*integral(f,aMin,bMin));

