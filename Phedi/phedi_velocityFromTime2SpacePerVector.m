function [t_mins_t_tip, x_mins_x_tip, PhotoLocation] = phedi_velocityFromTime2SpacePerVector(timeVec_4vel,RowOverTimeStruct,t_tips,Cf_PhediRegion)
% [t_mins_t_tip, x_mins_x_tip] = phedi_velocityFromTime2Space(timeVec_4vel,RowOverTimeStruct,t_tips)
% 
% OUTPUTS:
% this function returns a vector of t-t_tip and x-x_tip corrosponding to
% the vector timeVec_4vel.
% PhotoLocation - the approximated location where the phedis frame is
% located (in meters along the block);
%
% INPUTS:
% timeVec_4vel - time vector for corrosponding to the velocity vector,
%                   units of second.
% RowOverTimeStruct - output of the function 'IDT_PlotRowOverTime'
% t_tips -  optional but recomended. vector of t_tip for each phedi (in
%           units of seconds!!)
%
%% finding t-t_tip
t_mins_t_tip = bsxfun(@minus,repmat(timeVec_4vel(:),1,length(t_tips)),t_tips(:)');

%% creating location vector
x_mins_x_tip = -t_mins_t_tip*Cf_PhediRegion;

%% find location photographed
%--- find closest time sampled:
scratchRegionMeters = [0.07 0.09];
relevantRegion = RowOverTimeStruct.x>=min(scratchRegionMeters)&RowOverTimeStruct.x<=max(scratchRegionMeters);
[~, I] = min(abs(RowOverTimeStruct.frontTime_interp(relevantRegion)-mean(t_tips)));
relevantLocations = RowOverTimeStruct.x(relevantRegion);
PhotoLocation = relevantLocations(I);

end