function [Cf, frontAppearenceFrame, X, Y] = phantom_findCf_byIntensity(RowOverTime,frames2Search,varargin)
%Cf = phantom_findCf(RowOverTime)
%
% phantom_findCf will find the Cf in units of pixles/frame

[smth, UsrCutOff] = setDefaults4function(varargin,4,nan);

%% find columns on which to follow
%--- get typical row
meanRow = mean(RowOverTime(1:100,:),1);
%--- smooth and flip to find scratches
smoothedRow = smooth(meanRow,smth);
% upDwnRow = -smoothedRow;
upDwnRow = smoothedRow;
%--- find scratches by using find peaks
[~,locs,~,p] = findpeaks(upDwnRow);
% %--- choose only center third where intensity is high
% L = length(upDwnRow);
% borders = round([L/3, 2*L/3]);
% relevantPeak = (locs>=borders(1) & locs<=borders(2));
%--- choose all peaks:
relevantPeak = logical(ones(size(p)));
%--- take prominance of relevant peaks (actually scratches)
relevantP = p(relevantPeak);
relevantLocs = locs(relevantPeak);
%--- create cutoff for actual peaks (this assumes that there is noise!!)
if isnan(UsrCutOff)
%     CutOff = mean([min(relevantP), max(relevantP)]);
    CutOff = max(smoothedRow)/4;
else
    CutOff = UsrCutOff;
end
%--- take only actual peaks not noise
actualLogical = relevantP>CutOff;
actualLocs = relevantLocs(actualLogical);

%%
framesNums = 1:size(RowOverTime,1);
nanFrames = ~ismember(framesNums,frames2Search);

ROTDiff = bsxfun(@minus,RowOverTime,RowOverTime(1,:));
%--- create groups showing shifted movement
ROTDiffSaturLogic = ROTDiff>=1000 | ROTDiff<=-1000;
ROTDiffSaturNan = ROTDiffSaturLogic;
ROTDiffSaturNan(nanFrames,:) = 0;
CC = bwlabel(ROTDiffSaturNan);
%--- iterate on groups to find earliest point in each group
earlyCol = zeros(max(CC(:)),1);
earlyRow = zeros(max(CC(:)),1);
for Ic = 1:max(CC(:))
    [Crow,Ccol] = find(CC==Ic);
    [~,minRowIdx] = min(Crow(:));
    earlyCol(Ic) = Ccol(minRowIdx);
    earlyRow(Ic) = Crow(minRowIdx);
end
X = earlyCol;
Y = earlyRow;




% negativeLogical = ROTDiff(end,:)<0;
% RotLower = ROTDiff(:,negativeLogical);
% ABSrot = abs(RotLower);
% ABSrot(ABSrot>100)=100;
% ABSrot(nanFrames,:) = nan;
% [~,maxIdx] = max(ABSrot);
% 
% 
% maxIdx_relevant = maxIdx;
% % maxIdx_relevant(maxIdx_relevant<min(frames2Search) | maxIdx_relevant>max(frames2Search))=nan;
% lengthVec = 1:size(RowOverTime,2);
% X = find(negativeLogical); %lengthVec(~isnan(maxIdx_relevant));
% Y = maxIdx_relevant(~isnan(maxIdx_relevant));

[slope, frontAppearenceFrame ,~,~] = fitLinearModelWithOutliers(X,Y);
Cf = 1/slope;
