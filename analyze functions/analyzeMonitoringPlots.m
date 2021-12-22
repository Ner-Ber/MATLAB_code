function analyzeMonitoringPlots(DataStruct,varargin)
%ANALYZEMONITORINGPLOTS will display plots of the data needed to be reviewd
%in case of messup in the parameters
%
% optional figures: {'BigPicROT','sg','PhedisFront'}

if isempty(varargin)
    Plots = {'BigPicROT','sg','PhedisFront'};
else
    Plots = varargin;
end

%% plot IDT RowOverTime with edges
if any(strcmpi(Plots,'BigPicROT'))
    blockLength=0.15;
    figure; axBig = axes;
    IDT_PlotRowOverTime(DataStruct.BigPicRotStruct);
    hold on;
    plot([0 0],DataStruct.BigPicRotStruct.t([1,end]),'--g','LineWidth',2.0);
    plot([1 1]*blockLength,DataStruct.BigPicRotStruct.t([1,end]),'--g','LineWidth',2.0);
    
    %-- add small pic ROT
    spatialVec = DataStruct.AsperityData.spatialVec;
    timeCount = DataStruct.AsperityData.timeCount;
    RowOverTime = DataStruct.AsperityData.RowOverTime;
    PhotoLocation = DataStruct.PhediData.PhotoLocation;        %  in meteers on block
    DX = PhotoLocation-mean(spatialVec);
    spatialVecM = spatialVec+DX;
    smallROT = imagesc([spatialVecM(1) spatialVecM(end)],[timeCount(1) timeCount(end)],RowOverTime./max(RowOverTime(:)));

    uistack(smallROT,'bottom');
    uistack(smallROT,'up');
    
end
%% plot sg on pattern and its Gamma
if any(strcmpi(Plots,'sg'))
    gamma_plotGammaFit(DataStruct)
end

%% plot front from small camera
if any(strcmpi(Plots,'PhedisFront'))
    phantomPlotFrontFound(DataStruct);
end