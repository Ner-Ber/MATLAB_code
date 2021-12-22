function anls_plot_mesh(exp,event,Tstart,Tend,normalize,baseFig)

if nargin<6
    baseFig=gcf;
end

set(0,'DefaultFigureWindowStyle','docked')
e=acq132_event_get_data(exp,event,Tstart,Tend,1,'Sxy','Syy','Sxx','x_sg');
f=calc_Fx_stress_splined(e);

figure(baseFig);
baseFig=baseFig+1;
my_mesh(f.x(3:2:end-2),f.t,f.Sxy(:,3:2:end-2)./f.Syy(:,3:2:end-2),normalize)
title('Sxy/syy');
figure(baseFig);
baseFig=baseFig+1;
my_mesh(f.x(3:2:end-2),f.t,f.dSxxdx(:,2:2:end-1)./f.Syy(:,3:2:end-2),normalize);
title('dSxx/syy');
figure(baseFig);
baseFig=baseFig+1;
my_mesh(f.x(:,3:2:end-2),f.t,f.Fx(:,2:2:end-1)./f.Syy(:,3:2:end-2),normalize);
title('Fx/syy');

figure(baseFig);
baseFig=baseFig+1;
my_mesh(f.x(3:2:end-2),f.t,f.Sxy(:,3:2:end-2),normalize)
title('Sxy');
figure(baseFig);
baseFig=baseFig+1;
my_mesh(f.x(3:2:end-2),f.t,f.Syy(:,3:2:end-2),normalize)
title('Syy');
figure(baseFig);
baseFig=baseFig+1;
my_mesh(f.x(3:2:end-2),f.t,f.dSxxdx(:,2:2:end-1),normalize);
title('dSxx');
figure(baseFig);
baseFig=baseFig+1;
my_mesh(f.x(:,3:2:end-2),f.t,f.Fx(:,2:2:end-1),normalize);
title('Fx');