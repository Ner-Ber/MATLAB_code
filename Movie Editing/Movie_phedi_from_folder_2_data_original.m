function [analyzeAsperityCell, sg_data_cell, CAM_meta_cell, analyzePhediCell] = Movie_phedi_from_folder_2_data(varargin)
% [analyzeAsperityCell, sg_data_cell, CAM_meta_cell, analyzePhediCell] =
%       Movie_phedi_from_folder_2_data(preFrames, postFrames, experNum, eventNums_vec, res)
%
% Movie_phedi_from_folder_2_data will read the images from the movie, then
% claculate phedis and follow them end output cell arrays containig the
% relevant data.
%
% INPUTS:
% experNum, eventNums_vec, lineNum, res,          numOOM, preRowsTime, postRowsTime, prePhediTime, postPhediTime
% 1,         3,           'all',   32/(100e-6),   3,      0.5*1e-3,   4e-3,          1.5e-4,       1.5e-4
% are the defaults

%% set dafaults:
[experNum, eventNums_vec, lineNum, res, numOOM,...
    preRowsTime, postRowsTime, prePhediTime, postPhediTime] =...
    setDefaults4function(varargin,...
    1, 3, 'all', 32/(100e-6), 3,...
    0.5*1e-3, 4e-3, 1.5e-4, 1.5e-4);
eventNums_vec = eventNums_vec(:)';

%%

analyzeAsperityCell = cell(max(eventNums_vec),1);
sg_data_cell = cell(max(eventNums_vec),1);
CAM_meta_cell = cell(max(eventNums_vec),1);

analyzePhediCell = cell(max(eventNums_vec),1);

directoryList=my_dir;
exp_dir = directoryList{experNum};

for eventNum = eventNums_vec(:)'
    entireLoopTimer = tic;
    fprintf(1,['\n\nevent = ',num2str(eventNum),'\n']);
    [RowOverTimeCell, vallyPairsCell,~, ~,...
        maxTrajectories, CM_Trajectories, skewnessCell, varCell,center_of_massCell, sg_data, ~,Cam_Meta]...
        = analyzeEnlargementAndStrain(experNum, eventNum, preRowsTime, postRowsTime,lineNum);
    
    %---General stuff
    relevanrRowOverTime = RowOverTimeCell{end};
    fps = Cam_Meta.FrameRate;        % frames/s
    postFrames = round(postRowsTime*fps);
    preFrames = round(preRowsTime*fps);
    %--- regulate number of post/pre frames
    if postFrames>Cam_Meta.PostIms postFrames=Cam_Meta.PostIms; end
    if preFrames>(Cam_Meta.NumIms-Cam_Meta.PostIms) postFrames=Cam_Meta.PostIms; end
    
    if strcmp(Cam_Meta.CameraName,'PhIlya')
        expDetails=expDetailsRead(exp_dir);
        triggerDelay=expDetails.triggerDelay;
        fid= fopen([exp_dir '\Ph\' num2str(2) 'Time.bin'],'r');
        timeCount=(fread(fid,inf,'double')+triggerDelay);
        fclose(fid);
        eventFrame = Cam_Meta.NumIms-(Cam_Meta.PostIms+round(triggerDelay*Cam_Meta.FrameRate));
        if (eventFrame-preFrames)<1;    fromFrame=1; else	fromFrame=eventFrame-preFrames;	end
        if postFrames>Cam_Meta.PostIms;	toFrame=Cam_Meta.NumIms;	else	toFrame=eventFrame+postFrames;	end
        timeCount = timeCount(fromFrame:toFrame);
    else
        timeCount = [-preFrames:postFrames]/fps;
    end
    spatialVec = (0:(size(RowOverTimeCell{1},2)-1))/res;
    
    
    %--- calculate front velocity in scratches region
    BigPicRowOverTimeStruct = IDT_PlotRowOverTime(exp_dir,eventNum,preRowsTime,0,postRowsTime,0,1:Cam_Meta.ImageHeight,[],0);
    
    scratchRegionMeters = [0.07 0.09];
    xq = linspace(scratchRegionMeters(1),scratchRegionMeters(2),50);
    scartchLocationsInPix = interp1(BigPicRowOverTimeStruct.x,1:length(BigPicRowOverTimeStruct.x),xq);
    scartchLocationsInPix = unique([scartchLocationsInPix,round(min(scartchLocationsInPix)):round(max(scartchLocationsInPix))]);
    frontTimeInterped = interp1(BigPicRowOverTimeStruct.frontStepsPix,BigPicRowOverTimeStruct.frontStepsTime,scartchLocationsInPix);
    MeanFronVel = mean(diff(BigPicRowOverTimeStruct.x))*mean(diff(scartchLocationsInPix)./diff(frontTimeInterped))./mean(diff(BigPicRowOverTimeStruct.t));
    sg_data.Cf = MeanFronVel;
    
    %--- save asperity data
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
        analyzeAsperityStruct.postFrames,...
        sg_data_cell{eventNum},...
        CAM_meta_cell{eventNum}]...
        = deal(...
        RowOverTimeCell{end},...
        vallyPairsCell{end},maxTrajectories{end}, CM_Trajectories{end},...
        skewnessCell{end}, varCell{end},center_of_massCell{end}, spatialVec, res,...
        timeCount,preFrames, postFrames, sg_data, Cam_Meta);
    
    analyzeAsperityCell{eventNum} = analyzeAsperityStruct;
    
    
    
    %% calculate phedis
    disp('analyzing phedis...');
    
    useShiftedMethod = 1;
    if useShiftedMethod
        %--- calculate coarse shift
        [coarseShiftsVector, shiftedRowOverTime] = Movie_findCoarseShiftsInRowOverTime(relevanrRowOverTime,center_of_massCell{end});
        relevanrRowOverTime = shiftedRowOverTime;
    end
    
    %--- define frame to calc phedis on:
    prePhediFrames = round(prePhediTime*Cam_Meta.FrameRate);
    postPhediFrames = round(postPhediTime*Cam_Meta.FrameRate);
    firstFrame4Phedi = preFrames-prePhediFrames;
    lastFrame4Phedi = preFrames+postPhediFrames;
    
    % phedislocationInPix = phedisDB(2);
    smoothParameter = 0.001;
    phedisInitialLocsInPix = Movie_phedi_findTrenches(relevanrRowOverTime(firstFrame4Phedi,:), smoothParameter);
    %--- eliminate phedis close to edge:
    safetyDist = 100;
    phedisInitialLocsInPix = phedisInitialLocsInPix(abs(phedisInitialLocsInPix-length(relevanrRowOverTime(firstFrame4Phedi,:)))>safetyDist);
    phedisInitialLocsInPix = phedisInitialLocsInPix(phedisInitialLocsInPix>safetyDist);
    phedisInitialLocsInPix = phedisInitialLocsInPix(2:(end-1));
    
    calShiftTimer = tic;
    if useShiftedMethod
        [frameCount, PhediLocationPixShifted, slopeIncline, StrechingFactorsMat] = Movie_phedi_calculate_continues_shift(...
            relevanrRowOverTime, firstFrame4Phedi, lastFrame4Phedi, phedisInitialLocsInPix,numOOM,1);
        relevantCoarseShifts = coarseShiftsVector(firstFrame4Phedi:lastFrame4Phedi);
        PhediLocationPix = bsxfun(@plus,PhediLocationPixShifted,relevantCoarseShifts)-min(coarseShiftsVector(coarseShiftsVector>0));
        
        
    else
        %           [frameCount, PhediLocationPix, slopeIncline, StrechingFactorsMat] = Movie_phedi_calculate_continues_shift(...
        %               relevanrRowOverTime, firstFrame4Phedi, lastFrame4Phedi, phedisInitialLocsInPix,numOOM,0);
        [frameCount, PhediLocationPix, slopeIncline, StrechingFactorsMat] = Movie_phedi_calculate_continues_shift_updating(...
            relevanrRowOverTime, firstFrame4Phedi, lastFrame4Phedi, phedisInitialLocsInPix,numOOM,0);
    end
    T_shiftTimer = toc(calShiftTimer);
    disp(['time for calculating phedi shifts = ',num2str(T_shiftTimer)]);
    
    %--- transform units:
    PhediLocation = (PhediLocationPix-1)/res;
    measuredPhedisFromPlot = (phedisInitialLocsInPix-1)/res;
    timeVec= timeCount(frameCount);
    
    %--- save Phedi data
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
        ]...
        = deal(...
        PhediLocation, timeVec, measuredPhedisFromPlot, slopeIncline, firstFrame4Phedi, lastFrame4Phedi, phedisInitialLocsInPix, PhediLocationPix, StrechingFactorsMat);
    analyzePhediCell{eventNum} = analyzePhediStruct;
    
    T_iteration = toc(entireLoopTimer);
    disp(['Elapsed time for this event iteration is ',num2str(T_iteration),' seconds']);
end

end