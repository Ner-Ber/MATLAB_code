function [t_mins_t_tip, x_mins_x_tip, PhotoLocation] = phedi_velocityFromTime2Space(timeVec_4vel,RowOverTimeStruct,t_tips)
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


scratchRegionMeters = [0.07 0.09];
relevantRegion = RowOverTimeStruct.x>=min(scratchRegionMeters)&RowOverTimeStruct.x<=max(scratchRegionMeters);
Cf_scratchRegion = mean(RowOverTimeStruct.frontVel_interp(relevantRegion))*RowOverTimeStruct.res*RowOverTimeStruct.fps;

%% finding t_tip
if nargin<3
    t_tips = mean(RowOverTimeStruct.frontTime_interp(relevantRegion)); % assuming the imaging of the phedis is in the middle of the scratch region
    t_mins_t_tip = timeVec_4vel(:)-t_tips/RowOverTimeStruct.fps;
else
    t_mins_t_tip = timeVec_4vel(:)-t_tips;
end

%% creating location vector
x_mins_x_tip = -t_mins_t_tip*Cf_scratchRegion;

%% find location photographed
%--- find closest time sampled:
[~, I] = min(abs(RowOverTimeStruct.frontTime_interp(relevantRegion)-t_tips));
relevantLocations = RowOverTimeStruct.x(relevantRegion);
PhotoLocation = relevantLocations(I);

end