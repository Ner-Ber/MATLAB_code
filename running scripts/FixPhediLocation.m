function ExperimentCell = FixPhediLocation(PhediStructCell)
%% prep
%--- define the photoLocation:
defPhotoLoc = 0.08;

ExperimentCell = {};
for i=1:30
    if isfield(PhediStructCell{i},'theWorkspace')
        continue
    end
    %% take data from structure
    DataStruct = PhediStructCell{i};
    
    [x0,v0] = intersections(DataStruct.BigPicRotStruct.frontVelLoc_interpM,DataStruct.BigPicRotStruct.frontVel_interpMperS,...
        [1 1]*defPhotoLoc,[-1e5,1e5]);
    CfPixPerFrame = v0*DataStruct.ExperimentData.res./DataStruct.CamMeta.FrameRate;
    
    Phedidata = DataStruct.PhediData;
    AsperityData = DataStruct.AsperityData;
    BigPicRowOverTimeStruct = DataStruct.BigPicRotStruct;
    ExperimentData = DataStruct.ExperimentData;
    
    
    firstFrame4Phedi = Phedidata.firstFrame4Phedi;
    lastFrame4Phedi = Phedidata.lastFrame4Phedi;
    phedisInitialLocsInPix = Phedidata.phedisInitialLocsInPix;
    PhediLocation = Phedidata.PhediLocation;
    timeVec_4vel = Phedidata.timeVec_4vel_original;
    timeVec = Phedidata.timeVec;
    
    RowOverTime = AsperityData.RowOverTime;
    timeCount = AsperityData.timeCount;
    
    scratchRegionMeters = ExperimentData.scratchRegionMeters;
    
    %% find new front properties
    [intersectionFrame, tFrontPassing_Frms, crossLocInPix] = phedi_find_t_tips_wCf(...
        RowOverTime,CfPixPerFrame,phedisInitialLocsInPix,PhediLocation,firstFrame4Phedi:lastFrame4Phedi);
    %%
%     [tFrontPassing_Frms, crossLocInPix, CfPixPerFrm, frontAppearenceFrame] = phedi_findFrontTimeFrontSpeed(RowOverTime,firstFrame4Phedi:lastFrame4Phedi,phedisInitialLocsInPix,PhediLocation);
    t_tips = interp1(1:size(RowOverTime,1),timeCount,tFrontPassing_Frms);
%     meanFrame = mean(tFrontPassing_Frms);
%     t_phediCrossing = find_y0_at_x0(meanFrame,timeCount,1:size(RowOverTime,1));
    %--- find the corrosponding location for the time found
%     frontTime = (BigPicRowOverTimeStruct.frontTime_interp)/BigPicRowOverTimeStruct.fps;
%     potentialLocs = find_x0_at_y0(t_phediCrossing, frontTime, BigPicRowOverTimeStruct.x,'all');
    %--- find location closest to center of scratches
%     [~,I] = min(abs(potentialLocs-mean(scratchRegionMeters)));
    scratchCrossLocation = defPhotoLoc;
%     scratchCrossLocation = potentialLocs(I);
    if isempty(scratchCrossLocation)
        continue
    end
    %--- find crack velocity in this region
    deltaX = 0.001; %region in which to integrate velocity (in meters)
    relevantRegionLogical = BigPicRowOverTimeStruct.frontVelLoc_interpM>(scratchCrossLocation-deltaX) & BigPicRowOverTimeStruct.frontVelLoc_interpM<(scratchCrossLocation+deltaX);
    Cf_MetPerSec = mean(BigPicRowOverTimeStruct.frontVel_interpMperS(relevantRegionLogical));
    CfCalcFrom = 'big pic t_tips location';
    
    %--- calculate the Cf from the finding of the phedis region in big pic
    search_t_tips_region(1) = scratchRegionMeters(1)-0.02;
    search_t_tips_region(2) = scratchRegionMeters(2)+0.02;
%     searchRegionLogical = (search_t_tips_region(1)<=BigPicRowOverTimeStruct.x & search_t_tips_region(2)>=BigPicRowOverTimeStruct.x);
%     searchTimeAxis = BigPicRowOverTimeStruct.frontTime_interp/BigPicRowOverTimeStruct.fps;
%     searchTimeAxis(~searchRegionLogical) = nan;
%     [~,PhotoLocationIdx] = min(abs(searchTimeAxis-mean(t_tips)));
%     PhotoLocation = BigPicRowOverTimeStruct.x(PhotoLocationIdx);

    % [t_mins_t_tip_4vel, x_mins_x_tip_4vel, PhotoLocation] = phedi_velocityFromTime2Space(timeVec_4vel,BigPicRowOverTimeStruct,t_tips);
    [t_mins_t_tip_4vel, x_mins_x_tip_4vel, PhotoLocationUnused] =...
        phedi_velocityFromTime2SpacePerVector(timeVec_4vel,BigPicRowOverTimeStruct,t_tips,Cf_MetPerSec);
    
    try
        Xc = phedi_calcCohesiveZoneByVel(PhediVelocity,x_mins_x_tip_4vel);
    catch
        Xc = nan;
    end
    
    
    disp('saving phedi data...');
    analyzePhediStruct = struct(...
        'PhediLocation',PhediLocation,...
        'timeVec',timeVec,...
        'firstFrame4Phedi',firstFrame4Phedi,...
        'lastFrame4Phedi',lastFrame4Phedi,...
        'phedisInitialLocsInPix',phedisInitialLocsInPix,...
        'timeVec_4vel',timeVec_4vel,...
        't_mins_t_tip_4vel',t_mins_t_tip_4vel,...
        'x_mins_x_tip_4vel',x_mins_x_tip_4vel,...
        't_tips',t_tips,...
        'Xc',Xc,...
        'PhotoLocation',scratchCrossLocation,...
        'Cf',Cf_MetPerSec,...
        'CfCalcFrom',CfCalcFrom...
        );
    
    for fn = fieldnames(analyzePhediStruct)'
        Phedidata.(fn{1}) = analyzePhediStruct.(fn{1});
    end
    
    %--- update dataStruct
    DataStructNew = DataStruct;
    DataStructNew.PhediData = Phedidata;
    %--- smooth phedi signals:
%     DataStructNew = phedi_addSmoothedFunctions2Struct(DataStructNew);
    ExperimentCell{i} = DataStructNew;
end

save('fixedPhediStructCell_Exp5.mat','ExperimentCell','-v7.3');