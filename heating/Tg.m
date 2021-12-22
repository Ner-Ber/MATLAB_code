%[t,z,T]=Tg(dc,h,ArOverAn)
%
% tipically 0.01<ArOverAn<0.015 @ 5MPa (Dieterich & Kilgore 1994)
%       in my system with 150mm block its typically 6MPa 
% taur is the residual shear strain, typically ~2.5MPa

function [t,z,T]=Tg(dc,h,ArOverAn,taur)

v=logspace(-2,1,50);

alpha=1.1e-7; %Thermal diffusivity [m^2/s]
rho=1200;
c=1500; %heat capacity [J/kg/C]

% dc=10e-6;%1e-6 - 3e-6

tau_pr=1.6; %1 - 2 
tau_p=(2.5+tau_pr/2)*1e6; %factor 2 for mean value
A=5/333;% 0.01 - 0.02
%tauc=tau_p/A; % 150e6 - 350e6
%tauc=220e6;
% tauc=300e6;
tauc = taur./ArOverAn;

q=tauc*v;%rate of work
t=dc./v;%asperity life time
%t=[1e-7 1e-6 1e-5 1e-4 1e-3];
% h=dc/10;%0.3e-6;%layer height
% h=0;%0.3e-6;%layer height
%z=linspace(1e-10,h,1e3);
z=h;

%T=(v/(pi*alpha)).^0.5/(rho*c)*(dc^0.5*tauc); %Equivelent to z=0, (Rice 2006)

T=zeros(length(t),length(z));

for k=1:length(t)
    for l=1:length(z)
        f=@(tt) q(k)/2/(rho*c)*(pi*alpha*(t(k)-tt)).^(-0.5).*exp(-z(l).^2./(4*alpha*(t(k)-tt)));%factor 2 because we have two solids
        T(k,l)=integral(f,0,t(k));
    end
end


%T_avg=mean(T,2)+300;
T_avg=T+300;

% figure;
Name = ['tauc= ' num2str(round(tauc*1e-6)) ' dc= ' num2str(dc*1e6) ' h= ' num2str(h*1e6),' A/A_n=',num2str(ArOverAn),' \tau_r=',num2str(taur) ];
semilogx(v,T_avg,'-','DisplayName',Name);
% my_legend_add(['tauc= ' num2str(round(tauc*1e-6)) ' dc= ' num2str(dc*1e6) ' h= ' num2str(h*1e6),' A/A_n=',num2str(ArOverAn),' \tau_r=',num2str(taur) ]);
xlabel('v (m/s)');
ylabel('Tc (K)');

%figure(10);
% b=0.016;
% Ar=0.735-T_avg'/T_avg(1)*b.*log(v);
% % vm=100;
% % A0=0.65;
% % Ar=A0*( 1+T_avg'/T_avg(1)/A0*b.*log(1+vm./v) );
% semilogx(v,Ar,'-');
% my_legend_add(['tauc= ' num2str(round(tauc*1e-6)) ' dc= ' num2str(dc*1e6) ' h= ' num2str(h*1e6) ]);

% figure(2);
% a_b=0.001;
% vm=0.1e-6;
% taur=0.505 - a_b*T_avg'/T_avg(1).*log(1+v/vm);
% semilogx(v,taur,'-');
% my_legend_add(['tauc= ' num2str(round(tauc*1e-6)) ' dc= ' num2str(dc*1e6) ' h= ' num2str(h*1e6) ' a-b=' num2str(a_b) ]);

%j=2;plot(z/t(j)^0.5,T(j,:)/t(j)^0.5);

end