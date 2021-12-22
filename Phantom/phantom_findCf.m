function [Cf, frontAppearenceFrame, X, Y] = phantom_findCf(RowOverTime,frames2Search)
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

[slope, frontAppearenceFrame ,~,~] = fitLinearModelWithOutliers(X,Y);
Cf = 1/slope;
