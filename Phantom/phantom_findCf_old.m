function [Cf, frontAppearenceFrame, X, Y] = phantom_findCf(RowOverTime,frames2Search)
%Cf = phantom_findCf(RowOverTime)
%
% phantom_findCf will find the Cf in units of pixles/frame

framesNums = 1:size(RowOverTime,1);
nanFrames = ~ismember(framesNums,frames2Search);

RowOverTimeDiff = bsxfun(@minus,RowOverTime,RowOverTime(1,:));

% RowOverTimeSatur = RowOverTimeDiff;
% saturationLimit=50;
% RowOverTimeSatur(RowOverTimeSatur<-saturationLimit) = -saturationLimit;
% RowOverTimeSatur(RowOverTimeSatur>saturationLimit) = saturationLimit;
% 
% RowOverTimeDeriv = conv2(RowOverTimeSatur,[-1 1]','same');
% 
% RowOverTimeDerivNan = RowOverTimeDeriv;
% RowOverTimeDerivNan([1:10,end-10:end],:) = nan;
% 
% [~,maxIdx] = max(RowOverTimeDerivNan,[],1,'omitnan');

RowOverTimeSatur = abs(RowOverTimeDiff);
saturationLimit=300;
RowOverTimeSatur(RowOverTimeSatur<-saturationLimit) = -saturationLimit;
RowOverTimeSatur(RowOverTimeSatur>saturationLimit) = saturationLimit;

RowOverTimeSatur(nanFrames,:) = nan;

[~,maxIdx] = max(RowOverTimeSatur,[],1,'omitnan');
maxIdx_relevant = maxIdx;
% maxIdx_relevant(maxIdx_relevant<min(frames2Search) | maxIdx_relevant>max(frames2Search))=nan;
lengthVec = 1:size(RowOverTime,2);
X = lengthVec(~isnan(maxIdx_relevant));
Y = maxIdx_relevant(~isnan(maxIdx_relevant));

[slope, frontAppearenceFrame ,~,~] = fitLinearModelWithOutliers(X,Y);
Cf = 1/slope;
