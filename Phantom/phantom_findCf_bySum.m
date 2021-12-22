function [Cf, frontAppearenceFrame, X, Y] = phantom_findCf_bySum(RowOverTime,frames2Search, varargin)
%Cf = phantom_findCf(RowOverTime)
%
% phantom_findCf will find the Cf in units of pixles/frame

%% defaults
[Nbins] = setDefaults4function(varargin,10);
%% devide row over time into bins
N_rot = size(RowOverTime,2);
N_rotMod = N_rot-mod(N_rot,Nbins);
reshapedROT = reshape(RowOverTime(:,1:N_rotMod),size(RowOverTime,1),[],Nbins);
SummedROT = squeeze(sum(reshapedROT,2));

%% smooth the summed signals
SummedROTsmoothed = nan(size(SummedROT));
for i=1:Nbins
    SummedROTsmoothed(:,i) = supsmu(1:size(RowOverTime,1),SummedROT(:,i));
%     SummedROTsmoothed(:,i) = supsmu(1:size(RowOverTime,1),SummedROT(:,i),'Span',0.1);
end

%% second derivative
% ker = [-1 1]';
ker = [1 -2 1]';
SummedROT2ndDeriv = conv2(SummedROTsmoothed,ker,'same');

SummedROT2ndDerivSm = nan(size(SummedROT2ndDeriv));
for i=1:Nbins
    SummedROT2ndDerivSm(:,i) = supsmu(1:size(RowOverTime,1),SummedROT2ndDeriv(:,i));
end




%% old
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
