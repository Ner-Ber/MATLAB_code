function RowOverTimeStruct = phantom_PlotRowOverTime_BigPic(exp_dir,eventNum,varargin)
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
    -1.5e-3,1e-6,1.5e-3,0,'all',0);

% set(0,'DefaultFigureWindowStyle','docked')
% exp_details=expDetailsRead(exp_dir);
%--- set parameters
timeBase = 1; % 1 means in second

%% read stuff
ims=phantomGetLines_BigPic(exp_dir,eventNum,pre_time,intervalTime,post_time,timeBase,smt,'all',lineNum);

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

%--- calculate resolution
blockEdgesCorrdinates = [174 1076]; %[0 size(ims.lines,2)];
exp_details = expDetailsRead(exp_dir);
xx = linspace(0,exp_details.UpperBlockLength*1e-3,abs(diff(blockEdgesCorrdinates)));
res = abs(xx(2)-xx(1)); % meters per pixel
%--- spacial vector
alongFrameFixed = (1:size(ims.lines,2))-blockEdgesCorrdinates(1);
x = alongFrameFixed*res;

imagesc([x(1) x(end)],[ims.t(1) ims.t(end)],NormalizedRowOverTime);
% imagesc([ims.x(1) ims.x(end)],[ims.t(1) ims.t(end)],ims.lines);
ax = gca;
ax.YDir = 'normal';
title([ims.Date,'   ',ims.Exp,'   event=',num2str(ims.EventNum)]);
xlabel('x [m]'); ylabel('t [s]');
colorbar;
caxis([0.3 1.1]);

if nargout>0
    RowOverTimeStruct = struct;
    RowOverTimeStruct.DataMat = ims.lines;
    RowOverTimeStruct.DataMatNorm = NormalizedRowOverTime;
    RowOverTimeStruct.x = x;
    RowOverTimeStruct.t = ims.t;
    RowOverTimeStruct.details = [ims.Date,'   ',ims.Exp,'   event=',num2str(ims.EventNum)];
    RowOverTimeStruct.xlabel = 'x [m]';
    RowOverTimeStruct.ylabel = 't [s]';
end

end