function [Cf, frontAppearenceFrame, X, Y] = phantom_findCf_withGroups(RowOverTime,frames2Search)
%Cf = phantom_findCf(RowOverTime)
%
% phantom_findCf will find the Cf in units of pixles/frame

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

% AllCoors = [X(:),Y(:)];
% AllCoorsSoreted = sortrows(AllCoors);
% [~,ia,~] = unique(AllCoorsSoreted(:,1),'first') ;
% AllCoorsSelected = AllCoorsSoreted(ia,:);
% hold on; plot(AllCoorsSelected(:,1),AllCoorsSelected(:,2),'co');
% AllCoorsSoreted = sortrows(AllCoors);
% AllCoorsSoretedByRow = sortrows(fliplr(AllCoorsSelected));
% AllCoorsSoretedByRow = fliplr(sortrows(fliplr(AllCoorsSelected)));
% hold on; plot(AllCoorsSoretedByRow(1:2,1),AllCoorsSoretedByRow(1:2,2),'ro');
% hold on; plot(AllCoorsSoretedByRow(1:2,1),AllCoorsSoretedByRow(1:2,2),'r*');
% topPoints = AllCoorsSoretedByRow(1:2,:);


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
