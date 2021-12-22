function [A,T]=calcRate_State2(Cf_vec)


global vSlip;
global vSlip_t;
global d;

phi_0=5;%[sec]
phi_x=10^-8;%[sec]
beta=0.02;

Ar=Cf_vec*0;
vm=Cf_vec*0;

for j=1:length(Cf_vec)
    Cf=Cf_vec(j);
    %-----cohesive zone model
    sol=CrackSolutionGeneralCohesive(Cf,3.5e-7);
    vSlip_t=sol.t;%[sec]
    vSlip=-2*sol.v*sol.Uxx;%[m/s];
    d=sol.dc;
    %d=sol.dc/11; %note the factor 11 - use for dietrich law

    % %-----make a step function prifile
    %vSlip(vSlip_t>0)=max(vSlip);
    
%           [~,index]=min(abs(sol.x--30e-3));
%           vSlip(vSlip_t>0)=vSlip(index);
    
    %       vSlip(vSlip_t<0)=0+1e-5;
    
    %----solve Ode
    opts = odeset('RelTol',1e-4,'AbsTol',1e-8);
    [T,phi] = ode45( @phi_evol,[vSlip_t(1) vSlip_t(end)],phi_0,opts);
    A=1+beta*log(1+phi/phi_x);
    %A=1+beta*log(phi);
    
    slip=cumtrapz(vSlip)*mean(diff(sol.t));
    slipT=interp1(sol.t,slip,T);
    
    x=-T*sol.v;
    [~,index]=min(abs(x--30));
    Ar(j)=A(index)/A(1);
    %vm(j)=max(vSlip);
    
      [~,index]=min(abs(sol.x--30e-3));
      vm(j)=vSlip(index);
    %-----plot the results
    
%     figure(10);
%     hold all;
%     plot(slipT/d,A/A(1),'-');
%     xlim([0 7])
%     my_legend_add(num2str(sol.v/sol.Cr));
%     %
    figure(11);
    hold all;
    plot(-T*sol.v*1e3,A/A(1),'-');
    %plot(-T*sol.v*1e3,(A/A(1)-Ar(j))/(1-Ar(j)),'-');
    my_legend_add(num2str(sol.v/sol.Cr));
    my_xlim([-60 5]);
%     %
%     figure(12);
%     hold all;
%     plot(-vSlip_t*sol.v*1e3,vSlip);
%     my_legend_add(num2str(sol.v/sol.Cr));
end

% figure(13);
% semilogx(vm,Ar,'.-');


%---step v solution
%d=1.63E-6;
% V_x=d/phi_x;
% V=0.06;
% t=slipV_t;
% A=1+beta.*log(1+V_x/V+(phi_0/phi_x-V_x/V)*exp(-V*(t-t(1))/d));
%----------
% beta=0.04;
% d=1.63E-6;
% v=90e-3;
% phi_0=d/v*2;%[sec]
% tau=d/v;
% t=0:0.01e-4:0.4e-3;
% phi=(phi_0-tau)*exp(-t/tau)+tau;
% A=1+beta*log(phi);
% figure(4);
% plot(t*1e3,phi,'.-');
%  figure(3);
%  hold all;
%  plot(t*1e3,A/A(1),'.-');

function dy = phi_evol(t,phi)
global d;
vSlipTmp=calc_vSlip(t);
%dy = 1-(phi*vSlipTmp/2/d).^2;

%----Dietrich law
%dy = 1-phi*vSlipTmp/d;

%----slip law
dy = -phi*vSlipTmp/d*log(phi*vSlipTmp/d);

function vSlipTmp=calc_vSlip(t)
global vSlip;
global vSlip_t;
vSlipTmp=csaps(vSlip_t,vSlip,1,t);
%vSlipTmp=vSlip*(tanh(t/1e-9)+1)/2;
%v=0.06*ones(1,length(t));




