function acq132_stream_plot_loads(exp_dir,from,to,baseFig)
%exp_dir is the full path, including file name.   example:  '19-58-52\cal\stream2.mat'
% jump is used to dilute the data 30 is good

if (nargin<6)
    baseFig=1;
end
path=exp_dir;
acq132_stream_load_mat(path,from,to,'N','F');

%--------------------plot
set(0,'DefaultFigureWindowStyle','docked')

figure(baseFig);
baseFig=baseFig+1;
plot(t,F./N,'.-');
title('\mu');

% figure(baseFig);
% baseFig=baseFig+1;
% plot(t,F,'.-');
% title('Fs');
% 
% figure(baseFig);
% baseFig=baseFig+1;
% plot(t,N,'.-');
% title('FN');
 








