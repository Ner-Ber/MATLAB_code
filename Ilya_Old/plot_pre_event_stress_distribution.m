function plot_pre_event_stress_distribution(exp_dir,e_num,baseFig)
% e_num may be avector

if nargin <3
    baseFig=gcf;
end

ind=20;

for j=1:length(e_num)
e=acq132_event_get_data(exp_dir,e_num(j),'start','end',15,'Sxy','Sxx','Syy','x_sg','F','N');
e.Syy=my_smooth(e.Syy,3);
e.Sxx=my_smooth(e.Sxx,3);

spl=calc_Fx_stress_splined(e);
% mu=F(ind)/N(ind);
mu=0.5;

set(0,'DefaultFigureWindowStyle','docked');

figure(baseFig);
%baseFig=baseFig+1;
hold all
plot(spl.x(2:end-1),spl.Syy(ind,2:end-1),'o-');
title('Syy ');
%hold all
%plot(spl.x(2:end-1),spl.Fx(ind,:),'.-')
%title('Syy , Fx');
%  figure(baseFig+1);
%  plot(e.x_sg,e.Syy(ind,:),'o-')
%  hold all
 %plot(e.x_sg,e.Sxy(ind,:)./e.Syy(ind,:),'o-')
% title('Syy & Sxy/mu');

 
figure(baseFig+1);%tau/sigma
hold all
plot(spl.x(2:end-1),spl.Fx(ind,:)./spl.Syy(ind,2:end-1),'.-')
plot(e.x_sg,e.Sxy(ind,:)./e.Syy(ind,:),'.-');
%hold all
%plot(e.x_sg,e.Sxy(ind,:)./e.Syy(ind,:),'o-')
title('Fx/Syy')
end

%legend(cellstr(num2str(e_num')))

% figure(baseFig);
% title(['event -' num2str(e_num)])
% hold off
% figure(baseFig-1);
% hold off
%figure(baseFig+2);
