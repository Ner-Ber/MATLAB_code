function [DataStruct] = Movie_phedi_from_folder_2_data(varargin)
% [DataStruct] = Movie_phedi_from_folder_2_data(varargin)
%
% Movie_phedi_from_folder_2_data will read the images from the movie, then
% claculate phedis and follow them end output cell arrays containig the
% relevant data.
%
%INPUTS = DAFAULTS:
%1. experNum = 1
%2. eventNums_vec = 3
%3. lineNum = 'all'
%4. res = by phantom_getResFromScrathces function
%5. numOOM = 3
%6. preRowsTime = 3e-3
%7. postRowsTime = 3e-3
%8. prePhediTime = 1.5e-4
%9. postPhediTime = 3e-4

try
    entireAnalysisTimer = tic;
    %% set dafaults:
    [experNum, eventNum, lineNum, UsrRes, numOOM,...
        preRowsTime, postRowsTime, prePhediTime, postPhediTime] =...
        setDefaults4function(varargin,...
        1, 3, 'all', nan, 3,...
        3e-3, 3e-3, 1.5e-4, 3e-4);
    %% data on sg, cam_meta, and following asperities
    [RowOverTimeCell, vallyPairsCell,~, ~,...
        maxTrajectories, CM_Trajectories, skewnessCell, varCell,center_of_massCell, sgDataStruct, ~,CamMetaStruct]...
        = analyzeEnlargementAndStrain(experNum, eventNum, preRowsTime, postRowsTime,lineNum);
    
    %% General stuff
    RowOverTime = RowOverTimeCell{end};
    %--- get resolution
    if isnan(UsrRes)
        res = phantom_getResFromScrathces(RowOverTime);
    else
        res = UsrRes;
    end
    %--- get experiment details
    directoryList=my_dir;
    exp_dir = directoryList{experNum};
    
    %--- save general information
    [~,ExpDate,~] = fileparts(cd);
    ExpHour = exp_dir;
    
    ExperimentData = struct;
    
    [ExperimentData.experNum,...
        ExperimentData.eventNum,...
        ExperimentData.lineNum,...
        ExperimentData.res,...
        ExperimentData.numOOM,...
        ExperimentData.preRowsTime,...
        ExperimentData.postRowsTime,...
        ExperimentData.prePhediTime,...
        ExperimentData.postPhediTime,...
        ExperimentData.ExpDate,...
        ExperimentData.ExpHour,...
        ] = deal(experNum, eventNum, lineNum, res, numOOM,...
        preRowsTime, postRowsTime, prePhediTime, postPhediTime,ExpDate, ExpHour);
    
    
    
    center_of_mass = center_of_massCell{end};
    fps = CamMetaStruct.FrameRate;        % frames/s
    postFrames = round(postRowsTime*fps);
    preFrames = round(preRowsTime*fps);
    %--- regulate number of post/pre frames
    if postFrames>CamMetaStruct.PostIms; postFrames=CamMetaStruct.PostIms; end
    if preFrames>(CamMetaStruct.NumIms-CamMetaStruct.PostIms); postFrames=CamMetaStruct.PostIms; end
    
    if strcmp(CamMetaStruct.CameraName,'PhIlya')
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
    
    
    %--- calculate front velocity in scratches region
    BigPicRowOverTimeStruct = IDT_PlotRowOverTime(exp_dir,eventNum,preRowsTime,0,postRowsTime,0,1:CamMetaStruct.ImageHeight,[],[],0);
    %--- add velocities and vectors for each sg:
    sgDataStruct = sg_createVectorsForSg(sgDataStruct,BigPicRowOverTimeStruct);
    
    %% save asperity data
    disp('saving asperity data...');
    analyzeAsperityStruct = struct;
    [analyzeAsperityStruct.RowOverTime,...
        analyzeAsperityStruct.vallyPairs,...
        analyzeAsperityStruct.maxTrajectories,...
        analyzeAsperityStruct.CM_Trajectories,...
        analyzeAsperityStruct.skewness,...
        analyzeAsperityStruct.var,...
        analyzeAsperityStruct.center_of_mass,...
        analyzeAsperityStruct.spatialVec,...
        analyzeAsperityStruct.CameraResolution,...
        analyzeAsperityStruct.timeCount,...
        analyzeAsperityStruct.preFrames,...
        analyzeAsperityStruct.postFrames]...
        = deal(...
        RowOverTime,...
        vallyPairsCell{end},maxTrajectories{end}, CM_Trajectories{end},...
        skewnessCell{end}, varCell{end},center_of_mass, spatialVec, res,...
        timeCount,preFrames, postFrames);
    
    %% calculate phedis
    disp('analyzing phedis...');
    
    %--- define frame to calc phedis on:
    prePhediFrames = round(prePhediTime*CamMetaStruct.FrameRate);
    postPhediFrames = round(postPhediTime*CamMetaStruct.FrameRate);
    firstFrame4Phedi = preFrames-prePhediFrames;
    lastFrame4Phedi = preFrames+postPhediFrames;
    if firstFrame4Phedi<1; firstFrame4Phedi=1; end;
    if lastFrame4Phedi>size(RowOverTime,1); lastFrame4Phedi=size(RowOverTime,1); end;
    
    smoothParameter = 0.001;
    phedisInitialLocsInPix = Movie_phedi_findTrenches(RowOverTime(firstFrame4Phedi,:), smoothParameter);
    
    %--- eliminate phedis close to edge:
    safetyDist = 30;
    phedisInitialLocsInPix = phedisInitialLocsInPix(abs(phedisInitialLocsInPix-length(RowOverTime(firstFrame4Phedi,:)))>safetyDist);
    phedisInitialLocsInPix = phedisInitialLocsInPix(phedisInitialLocsInPix>safetyDist);
    phedisInitialLocsInPix = phedisInitialLocsInPix(2:(end-1));
    
    %--- eliminate close phedis
    scratchWidth = 100*1e-6;
    scratchWidth_inPix = scratchWidth*res;
    PhedisDist_inPix = abs(diff(phedisInitialLocsInPix));
    closePairs_logical = PhedisDist_inPix<0.35*scratchWidth_inPix;
    startings = find(([0;closePairs_logical]-[closePairs_logical;0])==-1);
    endings = find(([0;closePairs_logical]-[closePairs_logical;0])==1);
    discard = [];
    for k=1:length(startings); discard = cat(2,discard,(startings(k)+1):endings(k)); end
    phedisInitialLocsInPix(discard) = [];
    
    
    useShiftedMethod = 1;
    if useShiftedMethod
        %--- calculate coarse shift
        [coarseShiftsVector, shiftedRowOverTime] = Movie_findCoarseShiftsInRowOverTime(RowOverTime,center_of_mass);
        RowOverTime = shiftedRowOverTime;
    end
    
    
    calShiftTimer = tic;
    if useShiftedMethod
        [frameCount, PhediLocationPixShifted, slopeIncline, StrechingFactorsMat] = Movie_phedi_calculate_continues_shift(...
            RowOverTime, firstFrame4Phedi, lastFrame4Phedi, phedisInitialLocsInPix,numOOM,1);
        relevantCoarseShifts = coarseShiftsVector(firstFrame4Phedi:lastFrame4Phedi);
        if nnz(coarseShiftsVector>0)
            PhediLocationPix = bsxfun(@plus,PhediLocationPixShifted,relevantCoarseShifts)-min(coarseShiftsVector(coarseShiftsVector>0));
        else
            PhediLocationPix = bsxfun(@plus,PhediLocationPixShifted,relevantCoarseShifts);
        end
        
        
    else
        %           [frameCount, PhediLocationPix, slopeIncline, StrechingFactorsMat] = Movie_phedi_calculate_continues_shift(...
        %               relevanrRowOverTime, firstFrame4Phedi, lastFrame4Phedi, phedisInitialLocsInPix,numOOM,0);
        [frameCount, PhediLocationPix, slopeIncline, StrechingFactorsMat] = Movie_phedi_calculate_continues_shift_updating(...
            RowOverTime, firstFrame4Phedi, lastFrame4Phedi, phedisInitialLocsInPix,numOOM,0);
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
    % t_tips = phedi_findFrontTimeFromLocation(PhediLocation,timeVec);
    % [t_tips,Cf_pixPerSec] = phedi_findFrontTimeFromCenterOfMass...
    %     (center_of_mass,timeVec,firstFrame4Phedi:lastFrame4Phedi,analyzeAsperityStruct.RowOverTime,phedisInitialLocsInPix);
    % %--- calc Cf
    % Cf_MetPerSec = Cf_pixPerSec/res;
    
    [tFrontPassing_Frms , CfPixPerFrm, frontAppearenceFrame]= phedi_findFrontTimeFrontSpeed(...
        analyzeAsperityStruct.RowOverTime,firstFrame4Phedi:lastFrame4Phedi,phedisInitialLocsInPix,PhediLocationPix);
    t_tips = interp1(1:size(analyzeAsperityStruct.RowOverTime,1),timeCount,tFrontPassing_Frms);
    Cf_MetPerSec = CfPixPerFrm*fps/res;
    CfCalcFrom = 'scratches intensity drop';
    %--- find Cf by t_tips:
    scratchRegionMeters = [0.07 0.09];
    
    %--- calculate the Cf from the finding of the phedis region in big pic
    search_t_tips_region(1) = scratchRegionMeters(1)-0.02;
    search_t_tips_region(2) = scratchRegionMeters(2)+0.02;
    searchRegionLogical = (search_t_tips_region(1)<=BigPicRowOverTimeStruct.x & search_t_tips_region(2)>=BigPicRowOverTimeStruct.x);
    searchTimeAxis = BigPicRowOverTimeStruct.frontTime_interp/BigPicRowOverTimeStruct.fps;
    searchTimeAxis(~searchRegionLogical) = nan;
    [~,PhotoLocationIdx] = min(abs(searchTimeAxis-mean(t_tips)));
    PhotoLocation = BigPicRowOverTimeStruct.x(PhotoLocationIdx);
    PhediMaxDistBigPicPix = peak2peak(phedisInitialLocsInPix/res);  % distance between extreme phedis in big picture pixles
    %--- region to find Cf in:
    region4Cf = [floor(PhotoLocationIdx-2*ceil(PhediMaxDistBigPicPix)), ceil(PhotoLocationIdx+2*ceil(PhediMaxDistBigPicPix))];
    region4Cf(region4Cf<1) = 1;
    region4Cf(region4Cf<length(BigPicRowOverTimeStruct.frontVel_interpMperS)) = length(BigPicRowOverTimeStruct.frontVel_interpMperS);
    Cf_MetPerSecFromBigPic = mean(BigPicRowOverTimeStruct.frontVel_interpMperS(region4Cf(1):region4Cf(2)));
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
    analyzePhediStruct = struct;
    [analyzePhediStruct.PhediLocation,...
        analyzePhediStruct.timeVec,...
        analyzePhediStruct.measuredPhedisFromPlot,...
        analyzePhediStruct.slopeIncline,...
        analyzePhediStruct.firstFrame4Phedi,...
        analyzePhediStruct.lastFrame4Phedi,...
        analyzePhediStruct.phedisInitialLocsInPix,...
        analyzePhediStruct.PhediLocationPix,...
        analyzePhediStruct.StrechingFactorsMat,...
        analyzePhediStruct.PhediVelocity,...
        analyzePhediStruct.timeVec_4vel,...
        analyzePhediStruct.t_mins_t_tip_4vel,...
        analyzePhediStruct.x_mins_x_tip_4vel,...
        analyzePhediStruct.t_tips,...
        analyzePhediStruct.Xc,...
        analyzePhediStruct.PhotoLocation,...
        analyzePhediStruct.Cf,...
        analyzePhediStruct.CfCalcFrom...
        ] = deal(...
        PhediLocation, timeVec, measuredPhedisFromPlot, slopeIncline,...
        firstFrame4Phedi, lastFrame4Phedi, phedisInitialLocsInPix,...
        PhediLocationPix, StrechingFactorsMat,PhediVelocity, timeVec_4vel,...
        t_mins_t_tip_4vel,x_mins_x_tip_4vel,t_tips,Xc,...
        PhotoLocation,Cf_MetPerSec,CfCalcFrom...
        );
    
    
    %% calculate cohesive zone model
    disp('Calculating Gamma...');
    %--- find the gamma to fit to
    relevantSGs = find(scratchRegionMeters(1)<sgDataStruct.x_sg*1e-3 & sgDataStruct.x_sg*1e-3<scratchRegionMeters(2));
    relevantSG = round(mean(relevantSGs));  % take only on SG
    solAtSG = gamma_findGammaFit(sgDataStruct,BigPicRowOverTimeStruct, relevantSG);
    solAtSG.SG_calc = relevantSG;
    disp(['Gamma=',num2str(solAtSG.Gamma)]);
    fprintf('\n');
    %--- calc the cohesive zone model
    disp('calculating cohesive zone model on interface...');
    % [~, ~, Cr, ~ , ~, ~, ~,~,~,~]=CrackSolutionMaterialProperties;
    Cr = 1255;
    v = BigPicRowOverTimeStruct.AvgCf_scratches/Cr;
    h = 1e-7;
    CohesiveModelStruct=CrackSolutionGeneralCohesive_giveParameters(v,h,solAtSG.Gamma);
    % CohesiveModelStruct=CrackSolutionGeneralCohesive(v,h);
    CohesiveModelStruct.h = h;
    CohesiveModelStruct.v = v;
    CohesiveModelStruct.Cr = Cr;
    
    %% save and end
    DataStruct_0 = struct;
    [DataStruct_0.AsperityData,...
        DataStruct_0.PhediData,...
        DataStruct_0.SgData,...
        DataStruct_0.CamMeta,...
        DataStruct_0.ExperimentData,...
        DataStruct_0.BigPicRotStruct,...
        DataStruct_0.CohesiveModelStruct,...
        DataStruct_0.solAtSG...
        ] = deal(analyzeAsperityStruct,analyzePhediStruct,...
        sgDataStruct, CamMetaStruct,ExperimentData,BigPicRowOverTimeStruct,...
        CohesiveModelStruct,solAtSG);
    %--- smooth the data
    DataStruct = phedi_addSmoothedFunctions2Struct(DataStruct_0);
    T_iteration = toc(entireAnalysisTimer);
    disp(['Elapsed time for this event iteration is ',num2str(T_iteration),' seconds']);
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