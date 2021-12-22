function vSlip = calc_vSlip(v,sol)
%v=Uxx*Cf, function get Uxx*cf on y>0 and extrapolates to y=0 based on
%solSR


Gamma=2.5;%1.2;

%For Sub-Rayleigh
%v3_5=sol{2}.Uxx.max.*sol{2}.v'*sol{2}.Cr*(Gamma/sol{1}.Gamma)^0.5;
%v_0=sol{1}.Uxx.max.*sol{1}.v'*sol{1}.Cr*(Gamma/sol{1}.Gamma)^0.5;

%For super-shear
v3_5=sol{2}.Uxx.max.*sol{2}.v'*(Gamma/sol{1}.Gamma)^0.5;
v_0=sol{1}.Uxx.max.*sol{1}.v'*(Gamma/sol{1}.Gamma)^0.5;

%for j=1:length(v) vSlip{j}=interp1(v3_5,2*v_0,v{j}); end
vSlip=interp1(v3_5,2*v_0,v); 

%% 

%-plot
% for j=1:length(v) semilogx(vSlip{j},A{j},'ko'); hold all; end
% legend(fig.DisplayName);legend off;
% %ylim([0.5 0.9]);
% xlim([1e-2 30]);
% a=get(gca,'Children');
% for j=1:length(a) set(a(j),'color', fig.color{length(a)-j+1}); end    