function acq132_event_plot_ch_094_ex(exp_dir,event,pre,post,smt,baseFig)
% all the data is subtructed with the initial values. 


if (nargin<6)
    baseFig=1;
end

acq132_event_get_data(exp_dir,event,pre,post,smt,'ch_094_ex');

%--------------------plot
set(0,'DefaultFigureWindowStyle','docked')
figure(baseFig);
baseFig=baseFig+1;
plot(t,ch_094_ex(:,(3:5)),'.-');
title(['event-' num2str(event) 'ch-19,ch-23,ch-25']);

t_trig





