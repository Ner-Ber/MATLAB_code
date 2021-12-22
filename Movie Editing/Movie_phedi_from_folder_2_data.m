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
    %% set dafaults and general:
    [experNum, eventNum, lineNum, UsrRes, numOOM,...
        preRowsTime, postRowsTime, prePhediTime, postPhediTime] =...
        setDefaults4function(varargin,...
        1, 3, 'all', nan, 3,...
        3e-3, 3e-3, 1.5e-4, 3e-4);
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
    RowOverTime = phantom_getRowOverTime(exp_dir, eventNum, CamMetaStruct, preRowsTime, postRowsTime, lineNum);
    %% get resolution
    if isnan(UsrRes)
        res = phantom_getResFromScrathces(RowOverTime);
    else
        res = UsrRes;
    end
    
    %% create Experiment Data structure
    [~,ExpDate,~] = fileparts(cd);
    ExpHour = exp_dir;
    
    ExperimentData = struct(...
        'experNum',experNum,...
        'eventNum',eventNum,...
        'lineNum',lineNum,...
        'res',res,...
        'numOOM',numOOM,...
        'preRowsTime',preRowsTime,...
        'postRowsTime',postRowsTime,...
        'prePhediTime',prePhediTime,...
        'postPhediTime',postPhediTime,...
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
    BigPicRowOverTimeStruct = IDT_PlotRowOverTime(exp_dir,eventNum,preRowsTime,0,postRowsTime,0,1:CamMetaStruct.ImageHeight,[],[],0);
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
    prePhediFrames = round(prePhediTime*CamMetaStruct.FrameRate);
    postPhediFrames = round(postPhediTime*CamMetaStruct.FrameRate);
    firstFrame4Phedi = preFrames-prePhediFrames;
    lastFrame4Phedi = preFrames+postPhediFrames;
    if firstFrame4Phedi<1; firstFrame4Phedi=1; end;
    if lastFrame4Phedi>size(RowOverTime,1); lastFrame4Phedi=size(RowOverTime,1); end;
    
    %--- get phedis initial locations
    smoothParameter = 0.001;
    safetyDist = 30;
    scratchWidth = 100*1e-6;
    FirstRowSig = RowOverTime(firstFrame4Phedi,:);
    phedisInitialLocsInPix = phedi_getPhediInitialLocs(FirstRowSig,smoothParameter, safetyDist, res, scratchWidth);
    
    %--- calculate coarse shift
    [coarseShiftsVector, shiftedRowOverTime] = Movie_findCoarseShiftsInRowOverTime(RowOverTime,center_of_mass);
    %--- follow phedis
    calShiftTimer = tic;
    [frameCount, PhediLocationPixShifted, slopeIncline, StrechingFactorsMat] = Movie_phedi_calculate_continues_shift(...
        shiftedRowOverTime, firstFrame4Phedi, lastFrame4Phedi, phedisInitialLocsInPix,numOOM,1);
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
    % t_tips = phedi_findFrontTimeFromLocation(PhediLocation,timeVec);
    % [t_tips,Cf_pixPerSec] = phedi_findFrontTimeFromCenterOfMass...
    %     (center_of_mass,timeVec,firstFrame4Phedi:lastFrame4Phedi,RowOverTime,phedisInitialLocsInPix);
    % %--- calc Cf
    % Cf_MetPerSec = Cf_pixPerSec/res;
    
    [tFrontPassing_Frms , CfPixPerFrm, frontAppearenceFrame]= phedi_findFrontTimeFrontSpeed(...
        RowOverTime,firstFrame4Phedi:lastFrame4Phedi,phedisInitialLocsInPix,PhediLocationPix);
    t_tips = interp1(1:size(RowOverTime,1),timeCount,tFrontPassing_Frms);
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
    relevantSGs = find(scratchRegionMeters(1)<sgDataStruct.x_sg*1e-3 & sgDataStruct.x_sg*1e-3<scratchRegionMeters(2));
    relevantSG = round(mean(relevantSGs));  % take only on SG
    solAtSG = gamma_findGammaFit(sgDataStruct,BigPicRowOverTimeStruct, relevantSG);
    solAtSG.SG_calc = relevantSG;
    disp(['Gamma=',num2str(solAtSG.Gamma)]);
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