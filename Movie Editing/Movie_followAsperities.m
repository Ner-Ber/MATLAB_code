function [vallyPairs,asper_num,maxTrajectories,CM_Trajectories,skewnessMat,varMat,center_of_massMat, MarkedMatCell]...
    = Movie_followAsperities(RowOverTime)
% [vallyPairs,asper_num,maxTrajectories,CM_Trajectories,skewnessMat,varMat,center_of_massMat,MarkedMatCell]...
%    = Movie_followAsperities(RowOverTime)

%% Define the Asperities to follow
[vallyPairs,asper_num] = Movie_define_asperities_for_scratches(RowOverTime);
%% Make the images ready for analizing
MinOrMax = 'min';       % max to follow trench, 'min' to follow ridge
smthParameter = 1;
[MarkedMatCell, ~,RowOverTimeNormed] =...
    Movie_follow_maxima_in_row(RowOverTime, vallyPairs, smthParameter, MinOrMax);

%% Analyze the picture to get the trajectories of desired data:
%--- delete asperities of insifficient information due to noise: 
relevantAper = [];  % added 25/5/20 when trying to use a subset of all rows in frame. increases total noise. 
for strtPntIdx = 1:length(MarkedMatCell);   relevantAper=cat(1,relevantAper,~isnan(round(mean(find(MarkedMatCell{strtPntIdx}(1,:)))))); end;    % added 25/5/20 when trying to use a subset of all rows in frame. increases total noise. 
MarkedMatCell = MarkedMatCell(find(relevantAper));  % added 25/5/20 when trying to use a subset of all rows in frame. increases total noise. 


[maxTrajectories, CM_Trajectories, skewnessMat, varMat,center_of_massMat] =...
    Movie_analyze_asperities(MarkedMatCell, RowOverTimeNormed);

end