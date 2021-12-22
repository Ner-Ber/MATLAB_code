function [RowOverTime, timeVec] = phantom_getRowOverTime_and_timeCount(expPath, eventNum,varargin)
% [RowOverTime, timeVec] = phantom_getRowOverTime_and_timeCount(expPath, eventNum,varargin)
% varargin: pre_time, post_time, lineNum
% defaults: 0.5*1e-3, 4e-3,'all'

[pre_time, post_time, lineNum] = setDefaults4function(varargin,0.5*1e-3, 4e-3,'all');

[startl,interval,endl,timeVec]=phantomReadTime(expPath, eventNum,-pre_time,'min',post_time,1);
Images=phantomReadIms(expPath, eventNum,startl,interval,endl,1,lineNum);
ImagesMeaned = mean(Images,1);
RowOverTime = permute(ImagesMeaned,[3 2 1]);