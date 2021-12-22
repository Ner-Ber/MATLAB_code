function [maxTrajectories, CM_Trajectories, skewnessCell, varCell,center_of_massCell] = Movie_analyze_asperities(MarkedMatCell, RowOverTimeNormCell)
%[maxTrajectories, CM_Trajectories, skewnessCell, varCell,center_of_massCell] = Movie_analyze_asperities(MarkedMatCell, RowOverTimeNormCell)

[NumRelevantFrames, ImageWidth] = size(RowOverTimeNormCell{1});
xCoorMatrix = repmat(1:ImageWidth,NumRelevantFrames,1);
RowNums = length(MarkedMatCell);
maxTrajectories = {};
skewnessCell = {};
varCell = {};
center_of_massCell = {};
for rowIdx = 1:RowNums
    RowOverTimeNormed = RowOverTimeNormCell{rowIdx};
    maxTrajectories{rowIdx} = {};
    skewnessCell{rowIdx} = zeros(NumRelevantFrames,length(MarkedMatCell{rowIdx}));
    varCell{rowIdx} = zeros(NumRelevantFrames,length(MarkedMatCell{rowIdx}));
    for strtPntIdx = 1:length(MarkedMatCell{rowIdx});
        startingPixel = round(mean(find(MarkedMatCell{rowIdx}{strtPntIdx}(1,:))));
        scaningRange = 2;
        
        ridgeLogical = Movie_followRidge(RowOverTimeNormed, startingPixel, scaningRange, 'max');
        
        %% calculate maximum path:
        [y,x] = ind2sub(size(ridgeLogical),find(ridgeLogical));
        maxTrajectories{rowIdx}{strtPntIdx} = [x,y];
        
        %% calculate mean path:
        %--- find number of BW segmentation:
        relevantSegmentBW = MarkedMatCell{rowIdx}{strtPntIdx};
        %--- zero all places outside relevants:
        RowOverTimeZeroed = RowOverTimeNormed;
        RowOverTimeZeroed(~relevantSegmentBW) = 0;
        %--- compute mean!
        Center_of_Mass = sum(xCoorMatrix.*RowOverTimeZeroed,2)./sum(RowOverTimeZeroed,2);
        CM_Trajectories{rowIdx}{strtPntIdx} = [Center_of_Mass,sort(y)];
        
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
        skewnessCell{rowIdx}(:,strtPntIdx) = skewVec;
        varCell{rowIdx}(:,strtPntIdx) = varVec ;
        center_of_massCell{rowIdx}(:,strtPntIdx) = center_of_massVec ;
        
        
    end
end




end