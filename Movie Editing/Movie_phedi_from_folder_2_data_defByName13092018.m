function [DataStruct] = Movie_phedi_from_folder_2_data_defByName(experNum, eventNum,varargin)
% [DataStruct] = MOVIE_PHEDI_FROM_FOLDER_2_DATA_DEFBYNAME(experNum, eventNum,varargin)
%
% Movie_phedi_from_folder_2_data will read the images from the movie, then
% claculate phedis and follow them end output cell arrays containig the
% relevant data.
%
% Insert the optional parameters as follows:
%       'parameterName1',parameterVal1,'parameterName2',parameterVal2,...
%
%PARAMETERS = DAFAULTS:
% 'experNum',             1
% 'eventNum',             3
% 'lineNum',              'all'
% 'UsrRes',               by phantom_getResFromScrathces function
% 'numOOM',               3
% 'phediRelative2first',  1
% 'preRowsTime',          3e-3
% 'postRowsTime',         3e-3
% 'prePhediTime',         1.5e-4
% 'postPhediTime',        3e-4
% 'smooth4Phedi',         0.001
% 'PhediSafetyDist',      30    (in Pix)
% 'QuickAndDirty',        0
% 'scratchWidth',         100*1e-6
% 'scratchRegionMeters',  [0.07 0.09]
% 'Cr',                   1255
% 'CohsvModl_h',          1e-7
    

try
    entireAnalysisTimer = tic;
    %% set dafaults and general:
    
    
    
    defCell = {...
        'experNum',             1;...
        'eventNum',             3;...
        'lineNum',              'all';...
        'BlurROT',              20;...
        'UsrRes',               nan;...
        'numOOM',               3;...
        'phediRelative2first',  1;...
        'preRowsTime',          3e-3;...
        'postRowsTime',         3e-3;...
        'prePhediTime',         1.5e-4;...
        'postPhediTime',        3e-4;...
        'smooth4Phedi',         0.001;...
        'PhediSafetyDist',      30;...
        'QuickAndDirty',        0;...
        'scratchWidth',         100*1e-6;...
        'scratchRegionMeters',  [0.07 0.09];...
        'defPhotoLoc',          0.08;...
        'Cr',                   1255;...
        'CohsvModl_h',          1e-7...
        };
    Param = setDefaults4function_byName(defCell,varargin);
    
    %--- experiment path
    exper = my_dir;
    exp_dir = exper{experNum};
    
    %--- create the data struct
    DataStruct = struct;
    
    %% Camera mate data
    CamMetaStruct = CameraMetaAllCams(exp_dir);
    %--- update data structure
    DataStruct.CamMeta = CamMetaStruct;
    
    %% Get row over time
    RowOverTime = phantom_getRowOverTime(exp_dir, eventNum, CamMetaStruct, Param.preRowsTime, Param.postRowsTime, Param.lineNum);
    %--- blur rows of ROT if needed
    if Param.BlurROT>1 && mod(Param.BlurROT,1)==0 && isnumeric('none')
        g = createGaussianAproxCostume(1,Param.BlurROT);
        RowOverTime = conv2(RowOverTime,g,'same');
    end
    %% get resolution
    if isnan(Param.UsrRes)
        res = phantom_getResFromScrathces(RowOverTime);
    else
        res = Param.UsrRes;
    end
    
    %% create Experiment Data structure
    [~,ExpDate,~] = fileparts(cd);
    ExpHour = exp_dir;
    
    ExperimentData = struct(...
        'experNum',experNum,...
        'eventNum',eventNum,...
        'lineNum',Param.lineNum,...
        'res',res,...
        'numOOM',Param.numOOM,...
        'preRowsTime',Param.preRowsTime,...
        'postRowsTime',Param.postRowsTime,...
        'prePhediTime',Param.prePhediTime,...
        'postPhediTime',Param.postPhediTime,...
        'scratchRegionMeters',Param.scratchRegionMeters,...
        'ExpDate',ExpDate,...
        'ExpHour',ExpHour...
        );
    %--- update the data struct
    DataStruct.ExperimentData = ExperimentData;
    
    
    %% follow and analyze asperities
    disp('Definning and following asperities...');
    [vallyPairs,~,maxTrajectories,CM_Trajectories,skewnessMat,varMat,center_of_mass]...
        = Movie_followAsperities(RowOverTime);
    
    %%
    fps = CamMetaStruct.FrameRate;        % frames/s
    postFrames = round(Param.postRowsTime*fps);
    preFrames = round(Param.preRowsTime*fps);
    %--- regulate number of post/pre frames
    if postFrames>CamMetaStruct.PostIms; postFrames=CamMetaStruct.PostIms; end
    if preFrames>(CamMetaStruct.NumIms-CamMetaStruct.PostIms); postFrames=CamMetaStruct.PostIms; end
    
    if strcmp(CamMetaStruct.CameraName,'PhBig')
        expDetails=expDetailsRead(exp_dir);
        triggerDelay=expDetails.triggerDelay;
        fid= fopen([exp_dir '\Ph\' num2str(2) 'Time.bin'],'r');
        timeCount=(fread(fid,inf,'double')+triggerDelay);
        fclose(fid);
        eventFrame = CamMetaStruct.NumIms-(CamMetaStruct.PostIms+round(triggerDelay*CamMetaStruct.FrameRate));
        if (eventFrame-preFrames)<1;    fromFrame=1; else	fromFrame=eventFrame-preFrames;	end
        if postFrames>CamMetaStruct.PostIms;	toFrame=CamMetaStruct.NumIms;	else	toFrame=eventFrame+postFrames;	end
        timeCount = timeCount(fromFrame:toFrame);
    else
        timeCount = [-preFrames:postFrames]/fps;
    end
    spatialVec = (0:(size(RowOverTime,2)-1))/res;
    
    %% save asperity data structure
    disp('saving asperity data...');
    analyzeAsperityStruct = struct(...
        'RowOverTime',RowOverTime,...
        'vallyPairs',vallyPairs,...
        'maxTrajectories',{maxTrajectories},...
        'CM_Trajectories',{CM_Trajectories},...
        'skewness',skewnessMat,...
        'var',varMat,...
        'center_of_mass',center_of_mass,...
        'spatialVec',spatialVec,...
        'CameraResolution',res,...
        'timeCount',timeCount,...
        'preFrames',preFrames,...
        'postFrames',postFrames...
        );
    
    %--- update the data struct
    DataStruct.AsperityData = analyzeAsperityStruct;
    
    %% get Big Picture data
    SubFolders = my_dir(exp_dir);
    if (sum(strcmpi(SubFolders,'Ph'))+sum(strcmpi(SubFolders,'PhBig')))==2
        BigPicRowOverTimeStruct = phantom_BigPicStruct(exp_dir, eventNum,Param.preRowsTime,1e-10,Param.postRowsTime,1,1,'all');
    else
        BigPicRowOverTimeStruct = IDT_PlotRowOverTime(exp_dir,eventNum,Param.preRowsTime,0,Param.postRowsTime,0,1:CamMetaStruct.ImageHeight,[],[],0);
    end
    %--- update the data struct
    DataStruct.BigPicRotStruct = BigPicRowOverTimeStruct;
    
    
    %% get data from sg
    disp('Reading strain gages data...');
    sgDataStruct=acq132_event_get_data...
        (exp_dir,eventNum,'start','end',1,'Uxx','Uyy','Uxy','x_sg','y_sg','F','N');
    %--- add velocities and vectors for each sg and update DataStruct:
    sgDataStruct = sg_createVectorsForSg(sgDataStruct,BigPicRowOverTimeStruct);
    DataStruct.SgData = sgDataStruct;
    
    %% calculate phedis
    disp('analyzing phedis...');
    
    %--- define frame to calc phedis on:
    prePhediFrames = round(Param.prePhediTime*CamMetaStruct.FrameRate);
    postPhediFrames = round(Param.postPhediTime*CamMetaStruct.FrameRate);
    firstFrame4Phedi = preFrames-prePhediFrames;
    lastFrame4Phedi = preFrames+postPhediFrames;
    if firstFrame4Phedi<1; firstFrame4Phedi=1; end;
    if lastFrame4Phedi>size(RowOverTime,1); lastFrame4Phedi=size(RowOverTime,1); end;
    
    %--- get phedis initial locations
    FirstRowSig = RowOverTime(firstFrame4Phedi,:);
    phedisInitialLocsInPix = phedi_getPhediInitialLocs(FirstRowSig,Param.smooth4Phedi, Param.PhediSafetyDist, res, Param.scratchWidth);
    
    %--- calculate coarse shift
    [coarseShiftsVector, shiftedRowOverTime] = Movie_findCoarseShiftsInRowOverTime(RowOverTime,center_of_mass);
    %--- follow phedis
    calShiftTimer = tic;
    [frameCount, PhediLocationPixShifted, slopeIncline, StrechingFactorsMat] = Movie_phedi_calculate_continues_shift(...
        shiftedRowOverTime, firstFrame4Phedi, lastFrame4Phedi, phedisInitialLocsInPix,Param.numOOM,Param.phediRelative2first,Param.QuickAndDirty);
    relevantCoarseShifts = coarseShiftsVector(firstFrame4Phedi:lastFrame4Phedi);
    if nnz(coarseShiftsVector>0)
        PhediLocationPix = bsxfun(@plus,PhediLocationPixShifted,relevantCoarseShifts)-min(coarseShiftsVector(coarseShiftsVector>0));
    else
        PhediLocationPix = bsxfun(@plus,PhediLocationPixShifted,relevantCoarseShifts);
    end
    
    
    T_shiftTimer = toc(calShiftTimer);
    disp(['time for calculating phedi shifts = ',num2str(T_shiftTimer)]);
    
    %--- transform units: (minus 1 since the meters scale starts from 0 and the
    %pixel scale starts from 1)
    PhediLocation = (PhediLocationPix-1)/res;
    measuredPhedisFromPlot = (phedisInitialLocsInPix-1)/res;
    timeVec= timeCount(frameCount);
    
    %% phedi velocity
    [PhediVelocity, timeVec_4vel] = phedi_calcVelocity(PhediLocation,timeVec);
    
    %% phedi cohesive zone measurments
    %--- find crack velocity in this region
    deltaX = 0.001; %region in which to integrate velocity (in meters)
    relevantRegionLogical = BigPicRowOverTimeStruct.frontVelLoc_interpM>(Param.defPhotoLoc-deltaX) & BigPicRowOverTimeStruct.frontVelLoc_interpM<(Param.defPhotoLoc+deltaX);
    Cf_MetPerSec = mean(BigPicRowOverTimeStruct.frontVel_interpMperS(relevantRegionLogical));
    
    %--- find velocity in the photographed region
%     [x0,v0] = intersections(BigPicRowOverTimeStruct.frontVelLoc_interpM,BigPicRowOverTimeStruct.frontVel_interpMperS,...
%         [1 1]*Param.defPhotoLoc,[-1e5,1e5]);
    CfPixPerFrame = Cf_MetPerSec*DataStruct.ExperimentData.res./DataStruct.CamMeta.FrameRate;
    [intersectionFrame, tFrontPassing_Frms, crossLocInPix] = phedi_find_t_tips_wCf(...
        RowOverTime,CfPixPerFrame,phedisInitialLocsInPix,PhediLocation,firstFrame4Phedi:lastFrame4Phedi);
    t_tips = interp1(1:size(RowOverTime,1),timeCount,tFrontPassing_Frms);
    CfCalcFrom = 'big pic photo lcoation';
    
    % t_tips = phedi_findFrontTimeFromLocation(PhediLocation,timeVec);
    % [t_tips,Cf_pixPerSec] = phedi_findFrontTimeFromCenterOfMass...
    %     (center_of_mass,timeVec,firstFrame4Phedi:lastFrame4Phedi,RowOverTime,phedisInitialLocsInPix);
    % %--- calc Cf
    % Cf_MetPerSec = Cf_pixPerSec/res;
    
    %--- find approximate time of front passing phedis
%     [tFrontPassing_Frms, crossLocInPix, CfPixPerFrm, frontAppearenceFrame] = phedi_findFrontTimeFrontSpeed(RowOverTime,firstFrame4Phedi:lastFrame4Phedi,phedisInitialLocsInPix,PhediLocation);
%     t_tips = interp1(1:size(RowOverTime,1),timeCount,tFrontPassing_Frms);
    meanFrame = mean(tFrontPassing_Frms);
    t_phediCrossing = find_y0_at_x0(meanFrame,timeCount,1:size(RowOverTime,1));
    %--- find the corrosponding location for the time found
    frontTime = (BigPicRowOverTimeStruct.frontTime_interp)/BigPicRowOverTimeStruct.fps;
    potentialLocs = find_x0_at_y0(t_phediCrossing, frontTime, BigPicRowOverTimeStruct.x,'all');
    %--- find location closest to center of scratches
    [~,I] = min(abs(potentialLocs-mean(Param.scratchRegionMeters)));
    scratchCrossLocation = potentialLocs(I);
    %--- find crack velocity in this region
    deltaX = 0.001; %region in which to integrate velocity (in meters)
    relevantRegionLogical = BigPicRowOverTimeStruct.frontVelLoc_interpM>(scratchCrossLocation-deltaX) & BigPicRowOverTimeStruct.frontVelLoc_interpM<(scratchCrossLocation+deltaX);
    Cf_MetPerSec = mean(BigPicRowOverTimeStruct.frontVel_interpMperS(relevantRegionLogical));
%     CfCalcFrom = 'big pic t_tips location';
    
%     %-- data from big picture
%     Cf_IDT_phedi_MperS = find_y0_at_x0(mean(Param.scratchRegionMeters),...
%         BigPicRowOverTimeStruct.frontVel_interpMperS,BigPicRowOverTimeStruct.frontVelLoc_interpM);
%     Cf_Phantm_phedi_PixPerFrame = Cf_IDT_phedi_MperS*res/CamMetaStruct.FrameRate;
%     timeCrossing = find_y0_at_x0(mean(Param.scratchRegionMeters),...
%         BigPicRowOverTimeStruct.frontTime_interp/BigPicRowOverTimeStruct.fps,...
%         BigPicRowOverTimeStruct.x);
%     frame_phantom = find_y0_at_x0(timeCrossing,1:size(RowOverTime,1), timeCount);
    
%     [tFrontPassing_Frms , CfPixPerFrm, frontAppearenceFrame]= phedi_findFrontTimeFrontSpeed(...
%         RowOverTime,firstFrame4Phedi:lastFrame4Phedi,phedisInitialLocsInPix,PhediLocationPix);
%     [tFrontPassing_Frms , CfPixPerFrm, frontAppearenceFrame]= phedi_findFrontTimeFrontSpeed_wConstrains(...
%         RowOverTime,firstFrame4Phedi:lastFrame4Phedi,phedisInitialLocsInPix,PhediLocationPix,frame_phantom,Cf_Phantm_phedi_PixPerFrame);
%     t_tips = interp1(1:size(RowOverTime,1),timeCount,tFrontPassing_Frms);
%     Cf_MetPerSec = CfPixPerFrm*fps/res;
%     CfCalcFrom = 'scratches intensity drop';
    %--- find Cf by t_tips:
    
    %--- calculate the Cf from the finding of the phedis region in big pic
    search_t_tips_region(1) = Param.scratchRegionMeters(1)-0.02;
    search_t_tips_region(2) = Param.scratchRegionMeters(2)+0.02;
    searchRegionLogical = (search_t_tips_region(1)<=BigPicRowOverTimeStruct.x & search_t_tips_region(2)>=BigPicRowOverTimeStruct.x);
    searchTimeAxis = BigPicRowOverTimeStruct.frontTime_interp/BigPicRowOverTimeStruct.fps;
    searchTimeAxis(~searchRegionLogical) = nan;
    [~,PhotoLocationIdx] = min(abs(searchTimeAxis-mean(t_tips)));
    PhotoLocation = BigPicRowOverTimeStruct.x(PhotoLocationIdx);
%     PhediMaxDistBigPicPix = peak2peak(phedisInitialLocsInPix/res);  % distance between extreme phedis in big picture pixles
%     %--- region to find Cf in:
%     region4Cf = [floor(PhotoLocationIdx-2*ceil(PhediMaxDistBigPicPix)), ceil(PhotoLocationIdx+2*ceil(PhediMaxDistBigPicPix))];
%     region4Cf(region4Cf<1) = 1;
%     region4Cf(region4Cf<length(BigPicRowOverTimeStruct.frontVel_interpMperS)) = length(BigPicRowOverTimeStruct.frontVel_interpMperS);
%     Cf_MetPerSecFromBigPic = mean(BigPicRowOverTimeStruct.frontVel_interpMperS(region4Cf(1):region4Cf(2)));
    if Cf_MetPerSec==inf
        Cf_MetPerSec = Cf_MetPerSecFromBigPic;
        CfCalcFrom = 'big pic t_tips location';
    end
    % [t_mins_t_tip_4vel, x_mins_x_tip_4vel, PhotoLocation] = phedi_velocityFromTime2Space(timeVec_4vel,BigPicRowOverTimeStruct,t_tips);
    [t_mins_t_tip_4vel, x_mins_x_tip_4vel, PhotoLocationUnused] =...
        phedi_velocityFromTime2SpacePerVector(timeVec_4vel,BigPicRowOverTimeStruct,t_tips,Cf_MetPerSec);
    
    try
        Xc = phedi_calcCohesiveZoneByVel(PhediVelocity,x_mins_x_tip_4vel);
    catch
        Xc = nan;
    end
    
    %% save Phedi data
    disp('saving phedi data...');
    analyzePhediStruct = struct(...
        'PhediLocation',PhediLocation,...
        'timeVec',timeVec,...
        'measuredPhedisFromPlot',measuredPhedisFromPlot,...
        'slopeIncline',slopeIncline,...
        'firstFrame4Phedi',firstFrame4Phedi,...
        'lastFrame4Phedi',lastFrame4Phedi,...
        'phedisInitialLocsInPix',phedisInitialLocsInPix,...
        'PhediLocationPix',PhediLocationPix,...
        'StrechingFactorsMat',StrechingFactorsMat,...
        'PhediVelocity',PhediVelocity,...
        'timeVec_4vel',timeVec_4vel,...
        't_mins_t_tip_4vel',t_mins_t_tip_4vel,...
        'x_mins_x_tip_4vel',x_mins_x_tip_4vel,...
        't_tips',t_tips,...
        'Xc',Xc,...
        'PhotoLocation',PhotoLocation,...
        'Cf',Cf_MetPerSec,...
        'CfCalcFrom',CfCalcFrom...
        );
    
    %--- update dataStruct
    DataStruct.PhediData = analyzePhediStruct;
    %--- smooth phedi signals:
    DataStruct = phedi_addSmoothedFunctions2Struct(DataStruct);
    
    %% calculate cohesive zone model
    %--- turn off warning for this section
    warning('off');
    
    
    disp('Calculating Gamma...');
    %--- find the gamma to fit to
    relevantSGs = find(Param.scratchRegionMeters(1)<sgDataStruct.x_sg*1e-3 & sgDataStruct.x_sg*1e-3<Param.scratchRegionMeters(2));
    relevantSG = round(mean(relevantSGs));  % take only on SG
    solAtSG = gamma_findGammaFit(sgDataStruct,BigPicRowOverTimeStruct, relevantSG);
    solAtSG.SG_calc = relevantSG;
    disp(['Gamma=',num2str(solAtSG.Gamma)]);
    %--- calc the cohesive zone model
    disp('calculating cohesive zone model on interface...');
    % [~, ~, Cr, ~ , ~, ~, ~,~,~,~]=CrackSolutionMaterialProperties;

    v = BigPicRowOverTimeStruct.AvgCf_scratches/Param.Cr;
    CohesiveModelStruct=CrackSolutionGeneralCohesive_giveParameters(v,Param.CohsvModl_h,solAtSG.Gamma);
    % CohesiveModelStruct=CrackSolutionGeneralCohesive(v,h);
    CohesiveModelStruct.h = Param.CohsvModl_h;
    CohesiveModelStruct.v = v;
    CohesiveModelStruct.Cr = Param.Cr;
    
    
    %--- update dataStruct
    DataStruct.CohesiveModelStruct= CohesiveModelStruct;
    DataStruct.solAtSG= solAtSG;
    
    %--- turn back on warning messages
    warning('on');
    
    %% take time
    T_iteration = toc(entireAnalysisTimer);
    disp(['Elapsed time for this event iteration is ',num2str(T_iteration),' seconds']);
    fprintf('\n');
    fprintf('\n');
    
catch
    %--- save the workspace in case of error
    w = whos;
    for a = 1:length(w)
        theWorkspace.(w(a).name) = eval(w(a).name);
    end
    DataStruct = struct;
    DataStruct.theWorkspace = theWorkspace;
    putvar('DataStruct');
    
end

end