function [MarkedMatCell, AsperityDefinitionTraj, RowOverTimeNormed] =...
    Movie_follow_maxima_in_row(RowOverTime, vallyPairs, varargin)
%[MarkedMatCell, AsperityDefinitionTraj, RowOverTimeNormed] =...
%   Movie_follow_maxima_in_row(RowOverTime, vallyPairs, SmoothingParameter, MinOrMax)
%
%Movie_follow_maxima_in_row will create a trajectory of the average
%location of hieght distrubution of and asperity.
%
% smoothingParameter - defines the averging over each row-from-frame in
% order to reduce noise. Default is 7.

%% set dafaults:
%--- parameters cell array contains (by this order) the arguments: RelevantFrames,DarkThreshold,SmoothingParameter
[SmoothingParameter, MinOrMax] = setDefaults4function(varargin,7,'min');

%% followpeaks
%--- try to normalize the row over time before looking for ridges:
RowOverTimeShifted = bsxfun(@minus,RowOverTime,min(RowOverTime,[],2));
RowOverTimeNormed = bsxfun(@rdivide,RowOverTimeShifted,max(RowOverTimeShifted,[],2));
RowOverTimeNormed(isnan(RowOverTimeNormed)) = 0;

%--- smooth each frame
RowOverTimeSmoothed = conv2(RowOverTimeNormed,ones(1,SmoothingParameter)/SmoothingParameter,'same');

%--- find beggining points using vallyPairsCell:
MarkedMatCell = cell(size(vallyPairs,1),1);
AsperityDefinition = zeros(size(RowOverTimeSmoothed));
for begPairIdx = 1:size(vallyPairs,1)
    valleyA = Movie_followRidge(RowOverTimeSmoothed, min(vallyPairs(begPairIdx,:)), 2, MinOrMax);    % follow valley on one side
    valleyB = Movie_followRidge(RowOverTimeSmoothed, max(vallyPairs(begPairIdx,:)), 2, MinOrMax);    % follow valley on second side
    ValleysCombined = valleyB | valleyA;
    
    %--- mark area in which to define asperity
    currentMarkedMat = markBetweenOnes(ValleysCombined);
    MarkedMatCell{begPairIdx} = currentMarkedMat;
    AsperityDefinition = AsperityDefinition|ValleysCombined;
end

[I, J] = ind2sub(size(AsperityDefinition),find(AsperityDefinition));
AsperityDefinitionTraj = [I, J];


end

