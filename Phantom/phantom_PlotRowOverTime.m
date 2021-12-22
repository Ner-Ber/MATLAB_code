function RowOverTimeStruct = phantom_PlotRowOverTime(exp_dir,eventNum,varargin)
%RowOverTimeStruct = IDT_PlotRowOverTime(exp_dir,eventNum,pre_time,intervalTime,post_time,smt,lineNum,plotNew)
%will plot the row over time of the IDT camera
%MANDATORY
%exp_dir - path of the experiment folder
%eventNum - event number
%OPTIONAL:
%Tstart - starting time for display (trigger=0)
%Tinterval - interval in display frames
%Tend - interval in display frames
%smt - smooth factor along the x axis (default=0)
%lineNum - lines to consider (default = 'all');
%plotNew - open new figure? (def= 0)

%% set defaults
[pre_time,intervalTime,post_time,smt,lineNum,plotNew] = setDefaults4function(varargin,...
    -0.5e-3,1e-6,0.5e-3,0,'all',0);

% set(0,'DefaultFigureWindowStyle','docked')
% exp_details=expDetailsRead(exp_dir);
%--- set parameters
timeBase = 1; % 1 means in second

%% read stuff
ims=phantomGetLines(exp_dir,eventNum,pre_time,intervalTime,post_time,timeBase,smt,'all',lineNum);

%% plot stuff
% firRowMat = repmat(ims.lines(1,:),size(ims.lines,1),1);
% NormalizedRowOverTime = ims.lines./firRowMat;
NormalizedRowOverTime = bsxfun(@rdivide,ims.lines,ims.lines(1,:));

if plotNew
    figure;
else
    fig = gcf;
    figure(fig);
end
imagesc([ims.x(1) ims.x(end)],[ims.t(1) ims.t(end)],ims.lines)
title([ims.Date,'   ',ims.Exp,'   event=',num2str(ims.EventNum)]);
xlabel('x [m]'); ylabel('t [s]');

if nargout>0
    RowOverTimeStruct = struct;
    RowOverTimeStruct.DataMat = ims.lines;
    RowOverTimeStruct.DataMatNorm = NormalizedRowOverTime;
    RowOverTimeStruct.x = ims.x;
    RowOverTimeStruct.t = ims.t;
    RowOverTimeStruct.details = [ims.Date,'   ',ims.Exp,'   event=',num2str(ims.EventNum)];
    RowOverTimeStruct.xlabel = 'x [m]';
    RowOverTimeStruct.ylabel = 't [s]';
end

end