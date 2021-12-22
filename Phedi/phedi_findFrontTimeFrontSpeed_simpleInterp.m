function [tFrontPassing , CfPixPerFrm, frontAppearenceFrame] = phedi_findFrontTimeFrontSpeed_simpleInterp(RowOverTime,relevantFrames4Phedi,phediInitialLocs,PhediLocationPix)
%
%
% phedi_findFrontTimeFrontSpeed_simpleInterp will find the frames where the
% front is passing the phedis 


spatialVec = DataStruct.AsperityData.spatialVec;
timeCount = DataStruct.AsperityData.timeCount;
frames2Search = DataStruct.PhediData.firstFrame4Phedi:DataStruct.PhediData.lastFrame4Phedi;

[Cf, frontAppearenceFrame, Xint, Yint] = phantom_findCf(RowOverTime,frames2Search);
frontLoc = spatialVec(Xint);
frontTime = timeCount(Yint);
Cf_metPerSec = Cf*DataStruct.CamMeta.FrameRate/DataStruct.ExperimentData.res;