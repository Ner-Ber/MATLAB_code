function [Cf, frontAppearenceFrame, X, Y] = phantom_findCf_withCfBoundaries(RowOverTime,frames2Search,Cf_boundaries,time_boundaries)
%Cf = phantom_findCf(RowOverTime)
%
% phantom_findCf will find the Cf in units of pixles/frame

framesNums = 1:size(RowOverTime,1);
nanFrames = ~ismember(framesNums,frames2Search);

ROTDiff = bsxfun(@minus,RowOverTime,RowOverTime(1,:));
negativeLogical = ROTDiff(end,:)<0;
RotLower = ROTDiff(:,negativeLogical);
ABSrot = abs(RotLower);
ABSrot(ABSrot>100)=100;
ABSrot(nanFrames,:) = nan;
[~,maxIdx] = max(ABSrot);


maxIdx_relevant = maxIdx;
% maxIdx_relevant(maxIdx_relevant<min(frames2Search) | maxIdx_relevant>max(frames2Search))=nan;
lengthVec = 1:size(RowOverTime,2);
X = find(negativeLogical); %lengthVec(~isnan(maxIdx_relevant));
Y = maxIdx_relevant(~isnan(maxIdx_relevant));

%--- fit a model with contsrains
minSlope = 1/max(Cf_boundaries);
maxSlope = 1/min(Cf_boundaries);
min_t = min(time_boundaries);
max_t = max(time_boundaries);

P0 = [mean(1./Cf_boundaries), mean(time_boundaries)];
ub = [maxSlope, max_t];
lb = [minSlope, min_t];
model = @(P,x) P(1)*(x-size(RowOverTime,2))+P(2);

opts = optimset('Display','off');
P = lsqcurvefit(model,P0,X,Y,lb,ub,opts);
[slope, frontAppearenceFrame] = deal(P(1),P(2));
% [slope, frontAppearenceFrame ,~,~] = fitLinearModelWithOutliers(X,Y);
Cf = 1/slope;
