%[t,z,T]=Tg(dc,h,ArOverAn)
%
% tipically 0.01<ArOverAn<0.015 @ 5MPa (Dieterich & Kilgore 1994)
%       in my system with 150mm block its typically 6MPa
% taur is the residual shear strain, typically ~2.5MPa
%
% according to Dieterich&Kilgore 1996 a typical roughness for Acrylic
% (maybe like PMMA) will be in these values:
% ARation_vec = [1/430 4/500 16/730 1/340 4/390 16/520 1/400 4/400 16/460]=
%            =  [0.0023256        0.008     0.021918    0.0029412     0.010256     0.030769       0.0025         0.01     0.034783

function [v,Cf_vec,z,T_avg]=Tg_velocityProfile(Xc,h,ArOverAn,taur)
    
    Cf_vec = 10:10:1250;
    alpha=1.1e-7; %Thermal diffusivity [m^2/s]
    rho=1200;
    c=1500; %heat capacity [J/kg/C]
    tauc = taur./ArOverAn;
    z=h;
    %% load front velocities dependencies
    %-- peak velocity based on data
    %--***NOTICE!!! this is Cf=f(Vpeak)***
    CfVsVpeak = empModel_CfVsVpeak();
    
    %-- residual velocity based data
    %--***NOTICE!!! this is Vpeak=f(Cf)***
    VresVsCf = empModel_VresVsCf();
    
    %-- yield tau data
    %--***NOTICE!!! this is tau1st=f(Cf)***
%     TaucVsCf = empModel_Tau1stVsCf();
    %--***NOTICE!!! this is tau2nd=f(Cf)***
%     TaucVsCf = empModel_Tau2ndKinkVsCf();

    
    %     Xc = 3.5e-3;
    Tc = Xc./Cf_vec;
    T = nan(size(Cf_vec));
    %% iterate on various Cfs
    %-- create toy model
    Vintrp = linspace(0,2,5e2);
    corospndngCf = ppval(CfVsVpeak,Vintrp);
    [~,ia,~] = unique(corospndngCf);
    unqCf = corospndngCf(ia);
    unqVintrp = Vintrp(ia);
    for k=1:length(Cf_vec)
        Vpeak = interp1(unqCf,unqVintrp,Cf_vec(k));
        Vres = ppval(VresVsCf,Cf_vec(k));
        
        VofT = @(T) heaviside(Tc(k)/2-T).*heaviside(T).*(T.*2.*Vpeak./Tc(k)) +...
            heaviside(T-Tc(k)/2).*heaviside(Tc(k)-T).*((2*(Vres-Vpeak)./Tc(k)).*(T-Tc(k)/2)+Vpeak) + ...
            heaviside(T-Tc(k)).*Vres;
        
        
%         f=@(tt) ((ppval(TaucVsCf,Cf_vec(k))./ArOverAn).*VofT(tt)./2./((rho*c)*(pi*alpha*(Tc(k)-tt)).^(-0.5))).*exp(-z).^2./(4*alpha*(Tc(k)-tt));%factor 2 because we have two solids
        f=@(tt) ((tauc).*VofT(tt)./2./((rho*c)*(pi*alpha*(Tc(k)-tt)).^(-0.5))).*exp(-z).^2./(4*alpha*(Tc(k)-tt));%factor 2 because we have two solids
        T(k)=integral(f,0,10*Tc(k));
        
    end
    
    %T_avg=mean(T,2)+300;
    T_avg=T+273;
    v = interp1(unqCf,unqVintrp,Cf_vec);
    
%     figure;
    Name = ['tauc= ' num2str(round(tauc*1e-6)) ' X_c= ' num2str(Xc*1e6) ' h= ' num2str(h*1e6),' A/A_n=',num2str(ArOverAn),' \tau_r=',num2str(taur) ];
%     semilogx(v,T_avg,'-','DisplayName',Name);

    plot(Cf_vec,T_avg,'-','DisplayName',Name);
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