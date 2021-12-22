% phedi_getUxWithoutBigPic()
%
% will get phedi displacement for measurements without BigPic measurement
% (e.g. old measurements on lower block).


function [UxStruct] = phedi_getUxWithoutBigPic(experNum, eventNum, Param)
    %% get ROT
    %--- for phantom
    [RowOverTime, timeCount] = phantom_getRowOverTime_and_timeCount(experNum, eventNum, Param.preRowsTime, Param.postRowsTime);
    %--- for photron
    %     Pre = 5e3;
    %     Post = 5e3;
    %     [Images,Cam_Meta] = photronReadImages(experNum, eventNum,Pre,1,Post,1);
    %     RowOverTime = mean(Images,1);
    %     RowOverTime = squeeze(RowOverTime);
    %     RowOverTime = RowOverTime';
    %     Frms = -Pre:Post;
    %     timeCount = Frms/Cam_Meta.FrameRate;
    
    %% analyze
    [vallyPairs,~,maxTrajectories,CM_Trajectories,skewnessMat,varMat,center_of_mass, MarkedMatCell]...
        = Movie_followAsperities(RowOverTime);
    RowOverTime_normalized = ROT_normalizeIntesByAsperity(RowOverTime,MarkedMatCell);
    
    CamMetaStruct = CameraMetaAllCams(experNum);
    
    %%
    [~,centerTimeFrms] = min(abs(timeCount));
    firstFrame4Phedi = centerTimeFrms - round(Param.prePhediTime*CamMetaStruct.FrameRate);
    lastFrame4Phedi = centerTimeFrms + round(Param.postPhediTime*CamMetaStruct.FrameRate);
    
    %--- regulate number of post/pre frames
    spatialVec = (0:(size(RowOverTime,2)-1))/Param.res;
    RowOverTimePreBlur = RowOverTime;
    if Param.BlurROT>1 && mod(Param.BlurROT,1)==0 && isnumeric(Param.BlurROT)
        g = createGaussianAproxCostume(1,Param.BlurROT);
        RowOverTime = conv2(RowOverTime,g,'same');
    end
    
    disp('saving asperity data...');
    analyzeAsperityStruct = struct(...
        'RowOverTime',RowOverTime,...
        'RowOverTimePreBlur', RowOverTimePreBlur,...
        'RowOverTime_normalized',RowOverTime_normalized,...
        'vallyPairs',vallyPairs,...
        'maxTrajectories',{maxTrajectories},...
        'CM_Trajectories',{CM_Trajectories},...
        'skewness',skewnessMat,...
        'var',varMat,...
        'center_of_mass',center_of_mass,...
        'spatialVec',spatialVec,...
        'CameraResolution',Param.res,...
        'timeCount',timeCount);
    %         'preFrames',preFrames,...
    %         'postFrames',postFrames...
    %         );
    
    
    if firstFrame4Phedi<1; firstFrame4Phedi=1; end;
    if lastFrame4Phedi>size(RowOverTime,1); lastFrame4Phedi=size(RowOverTime,1); end;
    FirstRowSig = RowOverTime_normalized(firstFrame4Phedi,:);
    phedisInitialLocsInPix = phedi_getPhediInitialLocs(FirstRowSig,Param.smooth4Phedi, Param.PhediSafetyDist, Param.res, Param.scratchWidth);
    if Param.coarseShift
        [~, shiftedRowOverTime] = Movie_findCoarseShiftsInRowOverTime(RowOverTime_normalized,center_of_mass);
    else
        shiftedRowOverTime = RowOverTime_normalized;
        coarseShiftsVector = zeros(size(center_of_mass,1),1);
    end
    [frameCount, PhediLocationPixShifted, slopeIncline, StrechingFactorsMat, SlopeRegions] = Movie_phedi_calculate_continues_shift_total(...
        shiftedRowOverTime, firstFrame4Phedi, lastFrame4Phedi, phedisInitialLocsInPix,...
        Param.followMethod,Param.numOOM, Param.slopePenalty, Param.expandLowBound, Param.expandHighBound, Param.phediRelative2first);
    
    
    SpatialVec = (0:(length(FirstRowSig)-1))./Param.res;
    TemporalVecPix = firstFrame4Phedi:lastFrame4Phedi;
    TemporalVecSec = timeCount(TemporalVecPix);
    
    UxStruct = struct(...
        'RowOverTime', RowOverTime,...
        'res', Param.res,...
        'timeCount', timeCount,...
        'spacialVec', SpatialVec,...
        'TemporalVecPix', TemporalVecPix,...
        'TemporalVecSec', TemporalVecSec,...
        'frameCount', frameCount,...
        'PhediLocationPixShifted', PhediLocationPixShifted,...
        'slopeIncline', slopeIncline,...
        'StrechingFactorsMat', StrechingFactorsMat,...
        'SlopeRegions', SlopeRegions,...
        'coarseShiftsVector', coarseShiftsVector...
        );
    
    if Param.DoPlot
        figure; imagesc(UxStruct.spacialVec([1,end]),UxStruct.timeCount([1,end]),UxStruct.RowOverTime);
        hold on;
        plot((UxStruct.PhediLocationPixShifted-1)/UxStruct.res,UxStruct.TemporalVecSec,'-.');
    end
    
end