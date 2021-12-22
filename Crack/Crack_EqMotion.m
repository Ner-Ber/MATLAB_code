function [sol]=Crack_EqMotion(x,Uxy_0,Gamma)
%The function returns eq.of motion :Crack length , l, for each v;
%[v]=Cr;
%[tau_0]=Pa, Sxy_0-Sxy_r; (vector shoul be column oriented)
%[x]=m, (vector shoul be column oriented), should start from zero

[Cd Cs Cr, ~ , ~ ,E mu ~, ~ ,~]=CrackSolutionMaterialProperties;
k=Cs/Cd;%Broberg p.330
tau_0=Uxy_0*2*mu;


% % % %-------------Bi-lateral Crack
 %v=linspace(0.001,0.997,100);
% v=v*Cr;
%  R=@(z) (Cs^-2 - 2*z.^2).^2 + 4*z.^2 .* (Cd^-2 -z.^2).^0.5.* (Cs^-2-z.^2).^0.5; %just after eq. 6.3.42
%
% %-------Mode I
% % A=CrackSolutionCalcA(v,1);
% % for j=1:length(v)
% %     I_integrand=@(z) R(1i*z) .* (v(j)^-2 + z.^2).^-1.5 .* (Cd^-2 + z.^2).^-0.5;
% %     I(j)= (Cs^2/v(j) * integral( I_integrand, 0, Inf ) )^-1;
% % end
% % kappa=-I.* R(1./v) ./ ( Cs^-2* v.^-1 .*(v.^-2 - Cd.^-2).^0.5 );%eq6.367-eq6.368 is wrong
%
% % %------Mode II
%  A=CrackSolutionCalcA(v,2);
%  for j=1:length(v)
%      I_integrand=@(z) R(1i*z) .* (v(j)^-2 + z.^2).^-1.5 .* (Cs^-2 + z.^2).^-0.5;
%      I(j)= (Cs^2/v(j) *integral( I_integrand, 0, inf ) )^-1; %Sometimes doesn't converge dut to inf (although it decays as 1/z^2).  use the following form:
%     %I(j)= (Cs^2/v(j) * integral( I_integrand, 0, 1/v(j)*1e6 ) )^-1;
%  end
%  kappa=-I.* R(1./v) ./ ( Cs^-2* v.^-1 .*(v.^-2 - Cs.^-2).^0.5 );%eq6.367-eq6.368 is wrong
%  %-----
% kappa=real(kappa);%Due to numerical error beta has an imaginary part
% g=kappa.^2.*A;
% l=1/pi*Gamma*mu*4*(1-k^2)/tau_0^2./g;
% l0=1/pi*Gamma*mu*4*(1-k^2)/tau_0^2;


%%%------Semi-infinite crack -calculating S_+(h) Freund eq.6.4.18
%--------------ModeII
SPlus=@(v) CrackSolutionCalcSPlus(1./v,1/Cd,1/Cs);
A=@(v) CrackSolutionCalcA(v,2);
kappa=@(v) 1./SPlus(v).*(1-v/Cr)./(1-v/Cs).^0.5;%Dyanamic intensity factor Eq. 6.4.42
kappa=@(v) (1-v/Cr)./(1-v/Cs).^0.5;%Approximation
g=@(v) kappa(v).^2.*A(v);
%g=@(v) 1-v/Cr ; %Approximation p.356

%-----Mode III

% kappa=@(v) (1-v/Cs).^0.5;%Dyanamic intensity factor Eq. 6.4.42
% A=@(v) (1-(v/Cs).^2).^-0.5;
% g=@(v) kappa(v).^2.*A(v);

%-----------homogeneous tau_0
% l0=pi/8*Gamma*mu*4*(1-k^2)/tau_0^2; semi infinite crack
%  l0=0.8/pi*Gamma*mu*4*(1-k^2)/tau_0^2; %Edge crack. see our PNAS
%  v=linspace(0.001,0.997,100);
%  v=v*Cr;
%  l=l0./g(v);
%  sol.slipV=10^-3*v.*A(v);
% sol.SPlus=SPlus;



% %----------Clac KII_s for general profile of tau_0 
%x=[0;x];
%tau_0=[tau_0(1);tau_0];

tau=@(s) interp1(x,tau_0,s);
l=linspace(0,x(end),100);
v=l*0;
KII_s=l*0;
 
for j=1:length(l)
    
   F=@(s) 1+0.3*(1-s.^(5/4));
    I=@(s) 2*(1/pi/l(j))^0.5*tau(s).*F(s/l(j))./(1-(s/l(j)).^2 ).^0.5;  %Edge crack Tada p.197
    %I=@(s) 2*(1/pi/l(j))^0.5*tau(s)./(1-(s/l(j)).^2).^0.5;  %bi lateral crack Tada p.140
    %I=@(s) 2^0.5*(1/pi/l(j))^0.5*tau(s)./(1-s/l(j)).^0.5; % Semi infinte crack Freund p.352
    
    KII_s(j)=integral(I,0,l(j));
end
sol.KII_s=KII_s;


%calculate equtaion of motion for both edge crack and Eshelby
if (length(Gamma)==1)
    Gamma=l*0+Gamma;
else
    %G=[G(1);G];
    Gamma=interp1(x,Gamma,l);
end
sol.Gamma=Gamma;
sol.Gs=sol.KII_s.^2/(4*mu*(1-k^2)); %%Static energy release rate 

for j=1:length(l)
    if (Gamma(j)*4*mu*(1-k^2)-1*KII_s(j)^2<0) %Only above critical length
        fTmp= @(v)  Gamma(j)*4*mu*(1-k^2)-g(v)*KII_s(j)^2;
        v(j)=fzero(fTmp,v(j-1)+0.1);
    else
        v(j)=0;
        l0=l(j); %this is not accurate as it depends on resolution of vector l
    end
end



%%--------
sol.v=v/Cr;
sol.l=l; %[m]
sol.l0=l0;
sol.g=g;
sol.kappa=kappa;
sol.Cr=Cr;
sol.Cs=Cs;
sol.Cd=Cd;
sol.Gamma=Gamma;
sol.tau_0=tau_0;
sol.E=E;
sol.mu=mu;
sol.k=k;
sol.Uxy0=tau(l)/mu/2;

% % %-------Mode I
%  A=CrackSolutionCalcA(v,1);
%  kappa=1./SPlus.*(1-v/Cr)./(1-v/Cd).^0.5;%Dyanamic intensity factor Eq. 6.4.42
%g=kappa.^2.*A(v);
% %%-----
