function [PhediUxx_inTime] = phedi_strainFromVelocity(PhediVelocity,RowOverTimeStruct)
% [PhediUxx_inTime] = phedi_strainFromVelocity(PhediVelocity,RowOverTimeStruct)
% simple function to get local Uxx on the interface.
% PhediVelocity - uotput of the function 'phedi_calcVelocity'
% RowOverTimeStruct - uotput of the function 'IDT_PlotRowOverTime'


scratchRegionMeters = [0.07 0.09];
relevantRegion = RowOverTimeStruct.x>=min(scratchRegionMeters)&RowOverTimeStruct.x<=max(scratchRegionMeters);
Cf_scratchRegion = mean(RowOverTimeStruct.frontVel_interp(relevantRegion));

PhediUxx_inTime = -PhediVelocity./Cf_scratchRegion;

end