function [maxTrajectories, CM_Trajectories, skewnessMat, varMat,center_of_massMat] = Movie_analyze_asperities(MarkedMatCell, RowOverTimeNormed)
%[maxTrajectories, CM_Trajectories, skewnessCell, varCell,center_of_massCell] = Movie_analyze_asperities(MarkedMatCell, RowOverTimeNormCell)

[NumRelevantFrames, ImageWidth] = size(RowOverTimeNormed);
xCoorMatrix = repmat(1:ImageWidth,NumRelevantFrames,1);

CM_Trajectories = cell(size(MarkedMatCell));
maxTrajectories = cell(size(MarkedMatCell));
center_of_massMat = zeros(NumRelevantFrames,length(MarkedMatCell));
skewnessMat = zeros(NumRelevantFrames,length(MarkedMatCell));
varMat = zeros(NumRelevantFrames,length(MarkedMatCell));
for strtPntIdx = 1:length(MarkedMatCell);
    startingPixel = round(mean(find(MarkedMatCell{strtPntIdx}(1,:))));
    scaningRange = 2;
    ridgeLogical = Movie_followRidge(RowOverTimeNormed, startingPixel, scaningRange, 'max');
    
    %% calculate maximum path:
    [y,x] = ind2sub(size(ridgeLogical),find(ridgeLogical));
    maxTrajectories{strtPntIdx} = [x,y];
    
    %% calculate mean path:
    %--- find number of BW segmentation:
    relevantSegmentBW = MarkedMatCell{strtPntIdx};
    %--- zero all places outside relevants:
    RowOverTimeZeroed = RowOverTimeNormed;
    RowOverTimeZeroed(~relevantSegmentBW) = 0;
    %--- compute mean!
    Center_of_Mass = sum(xCoorMatrix.*RowOverTimeZeroed,2)./sum(RowOverTimeZeroed,2);
    CM_Trajectories{strtPntIdx} = [Center_of_Mass,sort(y)];
    
    %% calculate skewness and variance
    skewVec = zeros(NumRelevantFrames,1);
    varVec = zeros(NumRelevantFrames,1);
    center_of_massVec = zeros(NumRelevantFrames,1);
    for ii = 1:NumRelevantFrames
        considerOnly = ~~RowOverTimeZeroed(ii,:);
        
        xx=find(considerOnly);
        [center_of_mass,var_mass,skew_mass]=find_moments(xx,RowOverTimeZeroed(ii,considerOnly));
        center_of_massVec(ii) = center_of_mass;
        skewVec(ii) = skew_mass;
        varVec(ii) = var_mass;
        %             skewVec(ii) = skewness(RowOverTimeZeroed(ii,considerOnly),1);
        %             varVec (ii) = var(RowOverTimeZeroed(ii,considerOnly),1);
    end
    skewnessMat(:,strtPntIdx) = skewVec;
    varMat(:,strtPntIdx) = varVec ;
    center_of_massMat(:,strtPntIdx) = center_of_massVec ;
    
    
end

end