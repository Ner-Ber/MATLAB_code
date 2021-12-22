function [analyzePhediCell] = Movie_phedi_from_AsperCell_2_Phedis(analyzeAsperityCell,CAM_meta_cell, varargin)
% Movie_phedi_from_folder_2_data will read the images from the movie, then
% claculate phedis and follow them end output cell arrays containig the
% relevant data

%% set defaults
[numOOM ] = setDefaults4function(varargin,3);

%%

analyzePhediCell = cell(size(analyzeAsperityCell));

for eventNum = 1:length(analyzeAsperityCell)
    if isempty(analyzeAsperityCell{eventNum})
        continue
    end
    entireLoopTimer = tic;
    %% General stuff
    relevanrRowOverTime = analyzeAsperityCell{eventNum}.RowOverTime;
    fps = CAM_meta_cell{eventNum}.FrameRate;        % frames/s
    timeCount = analyzeAsperityCell{eventNum}.timeCount;
    res = analyzeAsperityCell{eventNum}.CameraResolution;
    spatialVec = (0:(size(relevanrRowOverTime,2)-1))/res;
    preFrames = analyzeAsperityCell{eventNum}.preFrames;
    postFrames = analyzeAsperityCell{eventNum}.postFrames;
    %% calculate phedis
    disp('analyzing phedis...');
    %--- define frame to calc phedis on:
    firstFrame4Phedi = preFrames-40;
    lastFrame4Phedi = preFrames+150;
    
    % phedislocationInPix = phedisDB(2);
    smoothParameter = 0.001;
    phedislocationInPix = Movie_phedi_findTrenches(relevanrRowOverTime(firstFrame4Phedi,:), smoothParameter);
    %--- eliminate phedis close to edge:
    safetyDist = 20;
    phedislocationInPix = phedislocationInPix(abs(phedislocationInPix-length(relevanrRowOverTime(firstFrame4Phedi,:)))>safetyDist);
    phedislocationInPix = phedislocationInPix(phedislocationInPix>safetyDist);
    phedislocationInPix = phedislocationInPix(2:(end-1));
    
    calShiftTimer = tic;
    [frameCount, PhediLocationPix, slopeIncline] = Movie_phedi_calculate_continues_shift(...
        relevanrRowOverTime, firstFrame4Phedi, lastFrame4Phedi, phedislocationInPix,...
        numOOM);
    T_shiftTimer = toc(calShiftTimer);
    disp(['tie for claculatin phedi shifts = ',num2str(T_shiftTimer)]);
    
    %--- transform units:
    PhediLocation = (PhediLocationPix-1)/res;
    measuredPhedisFromPlot = (phedislocationInPix-1)/res;
    timeVec= timeCount(frameCount);
    
    %--- save Phedi data
    disp('saving phedi data...');
    analyzePhediStruct = struct;
    [analyzePhediStruct.PhediLocation,...
        analyzePhediStruct.timeVec,...
        analyzePhediStruct.measuredPhedisFromPlot,...
        analyzePhediStruct.slopeIncline]...
        = deal(...
        PhediLocation, timeVec, measuredPhedisFromPlot, slopeIncline);
    analyzePhediCell{eventNum} = analyzePhediStruct;
    
    T_iteration = toc(entireLoopTimer);
    disp(['Elapsed time for this event iteration is ',num2str(T_iteration),' seconds']);
end

end