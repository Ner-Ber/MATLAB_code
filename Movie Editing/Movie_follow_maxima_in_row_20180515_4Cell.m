function [RowOverTimeNormCell,RowOverTimeCropedCell, RowOverTimeCell, MarkedMatCell, AsperityDefinitionTrajCell] =...
    Movie_follow_maxima_in_row(Images, vallyPairsCell, varargin)
%maximaTrajectories = Movie_follow_maxima_in_row(Images, vallyPairsCell, rowNumbers, RelevantFrames, DarkThreshold, SmoothingParameter, NumOfMaximas)
%
%Movie_follow_maxima_in_row will create a trajectory of the average
%location of hieght distrubution of and asperity.
%
% Images - a 3D matrix containing all frames of the movie in which you'll
% follow asperity.
% rowNums - numbers of rows to follow asperities in them.
% varargin will contain (by this order): rowNumbers, relevantFrames, darkThreshold,
% smoothingParameter, NumOfMaximas.
% 
% vallyPairsCell - a cell array, where each cell contains an nx2 matrix
% containing vally pairs, between then the asperity is defined. The
% function will follow the vadis, by this definning the asperity.
% numel(vallyPairsCell)==numel(rowNumbers) must hold!! if row numbers is
% not inserted vallyPairsCell=vallyPairsCell{1} will be.
%
% rowNumbers - row numbers to calculate. default is 30.
%
% relevantFrames - frames out of Images to cinsider in calculation. Default
% is all frames. if inputed as 'all', all frames will be considered.
%
% darkThreshold - a threshold under which all values will be considered as
% 0's. Dafault is 0.1/1.
%
% smoothingParameter - defines the averging over each row-from-frame in
% order to reduce noise. Default is 7.
%
% NumOfMaximas - an integer indicating the number of hieghest maximas to
% follow and calculate trajectoris of.  Default is 7.

%% set dafaults:
%--- parameters cell array contains (by this order) the arguments: RelevantFrames,DarkThreshold,SmoothingParameter
[DefRowNumbers, DefRelevantFrames, DefDarkThreshold, DefSmoothingParameter, DefNumOfMaximas, DefMinOrMax] =...
    deal(30, 'all', 0.1, 7, 7,'min');
Parameters = {DefRowNumbers, DefRelevantFrames, DefDarkThreshold, DefSmoothingParameter, DefNumOfMaximas, DefMinOrMax};
for i = 1:length(varargin)
    if ~isempty(varargin{i})
        Parameters{i} = varargin{i};
    end
end
[RowNumbers, RelevantFrames, DarkThreshold, SmoothingParameter, NumOfMaximas, MinOrMax] =...
    deal(Parameters{1},Parameters{2},Parameters{3},Parameters{4}, Parameters{5},Parameters{6});

vallyPairsCell = vallyPairsCell(1:min(length(RowNumbers),length(vallyPairsCell)));
RowNumbers = RowNumbers(1:min(length(RowNumbers),length(vallyPairsCell)));


[ImageHeight, ImageWidth, totNumOfFrames] = size(Images);

%% take relevant frames
if ischar(RelevantFrames)
    RelevantFrames = 1:totNumOfFrames;
end

%% followpeaks
%--- create cell array for each row:
MarkedMatCell = {};
RowOverTimeCell = {};
RowOverTimeNormCell = {};
RowOverTimeCropedCell = {};
%--- iterate over requested rows:
for rowIdx = 1:length(RowNumbers)
    begPair = vallyPairsCell{rowIdx};
    RowOverTime = Images(rowIdx,:,:);
    RowOverTime = permute(RowOverTime,[3 2 1]);
    
    %--- exclude complete-zero frames:
    rowSum = sum(RowOverTime,2);
    RelevantFrames = setdiff(RelevantFrames,find(~rowSum));
    
    %--- try to normalize the row over time before looking for ridges:
    RowOverTimeShifted = RowOverTime - repmat(min(RowOverTime,[],2),1,size(RowOverTime,2));
    RowOverTimeNormed = RowOverTimeShifted./repmat(max(RowOverTimeShifted,[],2),1,size(RowOverTime,2));
    RowOverTimeNormed(isnan(RowOverTimeNormed)) = 0;
    RowOverTimeSmoothed = RowOverTimeNormed;
    
    for frameIdx = RelevantFrames
        RowOverTimeSmoothed(frameIdx,:) = smooth(RowOverTimeSmoothed(frameIdx,:),SmoothingParameter);
    end
    
    RowOverTimeCroped=RowOverTimeSmoothed(RelevantFrames,:);
    
    
    %--- find beggining point using findpeaks:
%     %--- find all peaks of the forst frame in this row:
%     [pks,locs,w,p] = findpeaks(smooth(RowOverTimeCut(1,:),5));
%     %--- find hieghest peaks:
%     [~, IdniciesOfHieghest] = sort(pks, 'descend');
%     %--- relevant locations:
%     BeginLocs = locs(IdniciesOfHieghest(1:NumOfMaximas));
    
    %--- find beggining points using vallyPairsCell:
    BeginLocs = zeros(size(begPair,1),1);
    MarkedMatCell{rowIdx} = cell(size(begPair,1),1);
    AsperityDefinition = zeros(size(RowOverTimeCroped));
    for begPairIdx = 1:size(begPair,1)
        BeginLocs(begPairIdx) =round(mean(begPair(begPairIdx)));
        valleyA = Movie_followRidge(RowOverTimeCroped, min(begPair(begPairIdx,:)), 2, MinOrMax);    % follow valley on one side
        valleyB = Movie_followRidge(RowOverTimeCroped, max(begPair(begPairIdx,:)), 2, MinOrMax);    % follow valley on second side
        ValleysCombined = valleyB | valleyA;
        
        %--- mark area in which to define asperity
        currentMarkedMat = markBetweenOnes(ValleysCombined);
        MarkedMatCell{rowIdx}{begPairIdx} = currentMarkedMat;
        AsperityDefinition = AsperityDefinition|ValleysCombined;
    end

    RowOverTimeCell{rowIdx} = RowOverTime;
    RowOverTimeNormCell{rowIdx} = RowOverTimeNormed;
    RowOverTimeCropedCell{rowIdx} = RowOverTimeCroped;
    [I, J] = ind2sub(size(AsperityDefinition),find(AsperityDefinition));
    AsperityDefinitionTrajCell{rowIdx} = [I, J];
    
end


end

