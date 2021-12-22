function DataStructUpdated = Movie_fix_t_tips(DataStruct)
% DataStructUpdated = Movie_fix_t_tips(DataStruct)


if ~isfield(DataStruct.ExperimentData,'scratchRegionMeters')
    scratchRegionMeters = [0.07 0.09];
else
    scratchRegionMeters = DataStruct.ExperimentData.scratchRegionMeters;
end

RowOverTime = DataStruct.AsperityData.RowOverTime;
timeCount = DataStruct.AsperityData.timeCount;
firstFrame4Phedi = DataStruct.PhediData.firstFrame4Phedi;
lastFrame4Phedi = DataStruct.PhediData.lastFrame4Phedi;
phedisInitialLocsInPix = DataStruct.PhediData.phedisInitialLocsInPix;
PhediLocationPix = DataStruct.PhediData.PhediLocationPix;
res = DataStruct.ExperimentData.res;
FrameRate = DataStruct.CamMeta.FrameRate;

Cf_IDT_phedi_MperS = find_y0_at_x0(mean(scratchRegionMeters),...
    DataStruct.BigPicRotStruct.frontVel_interpMperS,DataStruct.BigPicRotStruct.frontVelLoc_interpM);
Cf_Phantm_phedi_PixPerFrame = Cf_IDT_phedi_MperS*res/FrameRate;
timeCrossing = find_y0_at_x0(mean(scratchRegionMeters),...
    DataStruct.BigPicRotStruct.frontTime_interp/DataStruct.BigPicRotStruct.fps,...
    DataStruct.BigPicRotStruct.x);
frame_phantom = find_y0_at_x0(timeCrossing,1:size(RowOverTime,1), timeCount);

[tFrontPassing_Frms , CfPixPerFrm, frontAppearenceFrame]= phedi_findFrontTimeFrontSpeed_wConstrains(...
    RowOverTime,firstFrame4Phedi:lastFrame4Phedi,phedisInitialLocsInPix,PhediLocationPix,frame_phantom,Cf_Phantm_phedi_PixPerFrame);
t_tips = interp1(1:size(RowOverTime,1),timeCount,tFrontPassing_Frms);
Cf_MetPerSec = CfPixPerFrm*FrameRate/res;

[t_mins_t_tip_4vel, x_mins_x_tip_4vel, PhotoLocationUnused] =...
        phedi_velocityFromTime2SpacePerVector(DataStruct.PhediData.timeVec_4vel,DataStruct.BigPicRotStruct,t_tips,Cf_MetPerSec);
    
%--- update the structure:
DataStructUpdated = DataStruct;
DataStructUpdated.PhediData.t_mins_t_tip_4vel = t_mins_t_tip_4vel;
DataStructUpdated.PhediData.x_mins_x_tip_4vel = x_mins_x_tip_4vel;
DataStructUpdated.PhediData.Cf = Cf_MetPerSec;
DataStructUpdated.PhediData.t_tips = t_tips;


end