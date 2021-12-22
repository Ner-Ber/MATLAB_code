function [N,F,FfromFx,NfromFx]=calcForcesFromStresses(acq)
%[N,F,FfromSxy]=calcForcesFromStresses(acq)

exp_dir=[acq.Path '\' acq.exp];
exp_details=expDetailsRead(exp_dir);
L=exp_details.UpperBlockLength;
W=exp_details.UpperBlockWidth;

spl.x=[0 acq.x_sg L];
%spl.Syy=spline(acq.x_sg,acq.Syy,spl.x);
spl.Syy=[acq.Syy(:,1) acq.Syy acq.Syy(:,end)];
N=trapz(spl.x,spl.Syy,2)*(1E6)*W*(1E-6)/9.8;%MPa->Pa,x,y coordinates are in mm , g

if(nargout>1)
%spl.Sxy=spline(acq.x_sg,acq.Sxy,spl.x);
spl.Sxy=[acq.Sxy(:,1) acq.Sxy acq.Sxy(:,end)];
F=trapz(spl.x,spl.Sxy,2)*(1E6)*W*(1E-6)/9.8;%MPa->Pa,x,y coordinates are in mm , g
end

if(nargout>2)
Fx=calc_Fx_stress_splined(acq);
spl.x=[0 Fx.x(2:end-1) L];
%spl.Fx=spline(Fx.x(2:end-1),Fx.Fx,spl.x);
spl.Fx=[Fx.Fx(:,1) Fx.Fx Fx.Fx(:,end)];
FfromFx=trapz(spl.x,spl.Fx,2)*(1E6)*W*(1E-6)/9.8;
end

if(nargout>3)
Fx=calc_Fx_stress_splined(acq);
spl.x=[0 Fx.x(2:end-1) L];
%spl.Fx=spline(Fx.x(2:end-1),Fx.Fx,spl.x);
spl.Syy0=[Fx.Syy0(:,1) Fx.Syy0 Fx.Syy0(:,end)];
NfromFx=trapz(spl.x,spl.Syy0,2)*(1E6)*W*(1E-6)/9.8;
end