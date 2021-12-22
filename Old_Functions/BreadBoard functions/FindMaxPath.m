function [pathsCoordinates,RidgesLogical] = FindMaxPath(IntensityMap,MinProm, MinWidth, GradThreshold)
%FindMaxPath will find and return the paths of local minimas of a row as
%their location moves down when propogating onward to lower (larger row
%number) y coordinates. 
%
%
% IntensityMap - is the m*n intensity map you want to follow maximas points
% on it.
%startingMaxPoints - is a vector containing starting points of local
%maximas to follow.

S = size(IntensityMap);

%% create matrices of peaks-per-row properties
% [PKS_pre,WIDTHS_pre,PROMS_pre] = deal(zeros(S));        % these are matrices containing data about peaks in their location according to the heat map
% for rowNum = 1:S(1)
%     [pks,locs,widths,proms] = findpeaks(IntensityMap(rowNum,:));
%     PKS_pre(rowNum,locs) = pks;
%     WIDTHS_pre(rowNum,locs) = widths;
%     PROMS_pre(rowNum,locs) = proms;
% end

%% select relevant peaks
% PromFilterLogical = PROMS_pre>=MinProm;
% WidthFilterLogical = WIDTHS_pre>=MinWidth;
% WP_filter = PromFilterLogical & WidthFilterLogical;     %combine width and prominence demands
% [PKS,WIDTHS,PROMS] = deal(PKS_pre,WIDTHS_pre,PROMS_pre);
% PKS(~WP_filter) = 0;
% WIDTHS(~WP_filter) = 0;
% PROMS(~WP_filter) = 0;
% 
% figure; imagesc(PKS);

% %% find traces of maximas
% Traces = zeros(S);
% for rowNum = 2:S(1)
%     
% end
% 

%% find initial peaks to begin with
[pks,locs,widths,proms] = findpeaks(IntensityMap(1,:));
%--- create filter for specifying the peaks
PromFilterLogical = proms>=MinProm;
WidthFilterLogical = widths>=MinWidth;
WP_filter = PromFilterLogical & WidthFilterLogical;     %combine width and prominence demands

%--- peaks after selection
PKS_selected = pks(WP_filter);
LOCS_selected = locs(WP_filter);

%% trace ridges on image
RidgesLogical = zeros(S);
pathsCoordinates = cell(length(LOCS_selected),1);
scaningRange = 2;
for PeakIdx = 1:length(LOCS_selected)
    startingPixel = LOCS_selected(PeakIdx);
    CurrentRidgeLogical = followRidge(IntensityMap, startingPixel, scaningRange, GradThreshold);
    
    %--- combine the ridges maps
    RidgesLogical = RidgesLogical|CurrentRidgeLogical;
    
    %--- save cooredinates of cuurent path
    [I, J] = ind2sub(size(CurrentRidgeLogical),find(CurrentRidgeLogical));
    pathsCoordinates{PeakIdx} = [I, J];
    
end





end