function phantomPlotFrontFound(DataStruct)

%% re-find the slope the front location
RowOverTime = DataStruct.AsperityData.RowOverTime;
spatialVec = DataStruct.AsperityData.spatialVec;
timeCount = DataStruct.AsperityData.timeCount;
frames2Search = DataStruct.PhediData.firstFrame4Phedi:DataStruct.PhediData.lastFrame4Phedi;

[Cf, frontAppearenceFrame, Xint, Yint] = phantom_findCf(RowOverTime,frames2Search);
X = spatialVec(Xint);
Y = timeCount(Yint);
Cf_metPerSec = Cf*DataStruct.CamMeta.FrameRate/DataStruct.ExperimentData.res;
% slope = 1/Cf;

%% plot
% frontFunc = @(x) frontAppearenceFrame+slope.*x;
appearTime = interp1(1:size(RowOverTime,1),timeCount,frontAppearenceFrame);
frontFunc = @(x) appearTime+(1./Cf_metPerSec).*x;

figure;
imagesc([spatialVec(1) spatialVec(end)],[timeCount(1) timeCount(end)],RowOverTime);
hold on;
plot(X,Y,'co');
% xx = 1:size(RowOverTime,2);
% plot(xx,frontFunc(xx),'r');
plot(spatialVec,frontFunc(spatialVec),'r');
title([DataStruct.BigPicRotStruct.details,'  Cf=',num2str(round(Cf_metPerSec)),'m/s']);

