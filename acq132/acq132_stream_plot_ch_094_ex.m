function acq132_stream_plot_ch_094_ex(exp_dir,from,to,baseFig)
%exp_dir is the full path, including file name.   example:  '19-58-52\cal\stream2.mat'
% all the data is subtructed with the initial values. 


if (nargin<4)
    baseFig=gcf;
end
path=exp_dir;
acq132_stream_load_mat(path,from,to,'ch_094_ex');


%---offset correction to ch_094_ex-notice it's done after the data cut
ch_094_ex(:,2:5)=subtruct_norm(ch_094_ex(:,2:5));

%--------------------plot
set(0,'DefaultFigureWindowStyle','docked')
figure(baseFig);
baseFig=baseFig+1;
plot(t,ch_094_ex(:,3:5),'.-');

title('ch-19,ch-23,ch-25');
my_xlim([from to])





