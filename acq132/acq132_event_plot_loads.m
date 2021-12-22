function acq132_event_plot_loads(exp_dir,event,pre,post,smt,baseFig)

sens4=18.37/1000;
sens5=11.66/1000;
sensF=2.94;%2.26;

if (nargin<5)
    baseFig=1;
end

acq132_event_get_data(exp_dir,event,pre,post,smt,'F','ch_094_ex');

FirstLine=repmat(mean(ch_094_ex(1:20,2:end)),length(ch_094_ex(:,1)),1);
ch_094_ex(:,2:end)=ch_094_ex(:,2:end)-FirstLine;
ch_094_ex(:,4)=ch_094_ex(:,4)/sens4;
ch_094_ex(:,5)=ch_094_ex(:,5)/sens5;

t_trig
%--------------------plot
set(0,'DefaultFigureWindowStyle','docked')

figure(baseFig);
baseFig=baseFig+1;
plot(t,ch_094_ex(:,2),'.-');
hold all
plot(t,ch_094_ex(:,3),'.-'); 
plot(t,ch_094_ex(:,3)-smooth(ch_094_ex(:,2),41),'.-');
plot(t,ch_094_ex(:,3)+(F-mean(F(1:20)))/sensF,'.-');
title('ch-17,ch-19');
hold off

figure(baseFig);
baseFig=baseFig+1;
plot(t,ch_094_ex(:,4:5),'.-');
hold all
plot(t,ch_094_ex(:,4:5)-repmat(smooth(ch_094_ex(:,2),41),1,2),'.-');
plot(t,ch_094_ex(:,4:5)+repmat((F-mean(F(1:20)))/sensF,1,2),'.-');
title('ch-23,ch-25');
hold off

figure(baseFig);
baseFig=baseFig+1;
plot(smooth(ch_094_ex(:,2),51),F-mean(F(1:20)),'.-');
title('Fs');








