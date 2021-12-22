function NewDataStruct = phedi_addOrChngPhediTtips(DataStruct)
    
    Param = DataStruct.Param;
    BigPicRotStruct = DataStruct.BigPicRotStruct;
    PhediLocation = DataStruct.PhediData.PhediLocation;
    measuredPhedisFromPlot = DataStruct.PhediData.measuredPhedisFromPlot;
    timeVec = DataStruct.PhediData.timeVec;
    
    %% get photo location and time
    expMeasures = expMeasureProperties(DataStruct.ExperimentData.ExpHour);
    PhotoLocation = expMeasures.photoLocation;
    
    %% find phedis location and time
    [x_tips,t_tips] = phedi_find_t_tips_wPhotoLocation(BigPicRotStruct,PhotoLocation,DataStruct.AsperityData.spatialVec,...
        PhediLocation,measuredPhedisFromPlot,timeVec);
    
    
    %% phedi t_tips
    %--- find crack velocity in this region using big picture
    deltaX = 0.001; %region in which to integrate velocity (in meters)
    relevantRegionLogical = BigPicRotStruct.frontVelLoc_interpM>(Param.defPhotoLoc-deltaX) & BigPicRotStruct.frontVelLoc_interpM<(Param.defPhotoLoc+deltaX);
    Cf_MetPerSec = mean(BigPicRotStruct.frontVel_interpMperS(relevantRegionLogical));
    
    %--- get t_tips
    %             CfPixPerFrame = Cf_MetPerSec*DataStruct.ExperimentData.res./DataStruct.CamMeta.FrameRate;
    %             [intersectionFrame, tFrontPassing_Frms, crossLocInPix] = phedi_find_t_tips_wCf(...
    %                 RowOverTime,CfPixPerFrame,phedisInitialLocsInPix,PhediLocation,firstFrame4Phedi:lastFrame4Phedi);
    %             t_tips = interp1(1:size(RowOverTime,1),timeCount,tFrontPassing_Frms);
    CfCalcFrom = 'big pic photo lcoation';
    
    
    %--- get t_mins_t_tip
    t_mins_t_tip = bsxfun(@minus,timeVec(:),t_tips(:)');
    %% save the new data
    NewDataStruct = DataStruct;
    NewDataStruct.PhediData.CfCalcFrom = CfCalcFrom;
    NewDataStruct.PhediData.t_mins_t_tip = t_mins_t_tip;
    NewDataStruct.PhediData.Cf = Cf_MetPerSec;
    NewDataStruct.PhediData.x_tips = x_tips;
    NewDataStruct.PhediData.t_tips = t_tips;
    
end