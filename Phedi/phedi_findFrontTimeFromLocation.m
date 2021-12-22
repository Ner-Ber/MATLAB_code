function tFrontPassing = phedi_findFrontTimeFromLocation(phediLocationMatrix,phediTimeVector,varargin)
% tFrontPassing = phedi_findFrontTimeFromLocation(phediLocationMatrix,phediTimeVector)
%
% phedi_findFrontTimeFromLocation - will approximate the time of the
% passing front by the asumption that there is no movement of the phedi
% before the front. Movement is detected by a threshold (NoiseFactor)
% comparing to the noise in the begining of the measurment.
%
% VARARGIN = DEFAULTS:
% timeBeforEvent = -11e-5	(time to sample noise)
% N_firstSteps = 40         (steps to sample noise in case time is irrelevant)
% NoiseFactor = 8;          (maximal noise time this will determine first movement)

%% set defaults:
[timeBeforEvent, N_firstSteps, NoiseFactor] = setDefaults4function(varargin,-11e-5,40,8);

phediLocationMatrix = bsxfun(@minus, phediLocationMatrix,phediLocationMatrix(1,:));
%% get noise of location before event
NoiseCrop = phediLocationMatrix(phediTimeVector<timeBeforEvent,:);
NoiseCrop_time = phediTimeVector(phediTimeVector<timeBeforEvent);
if isempty(NoiseCrop_time)
    NoiseCrop = phediLocationMatrix(1:N_firstSteps,:);
    NoiseCrop_time = phediTimeVector(1:N_firstSteps);
end

maximalNoise = max(abs(NoiseCrop),[],1);

%% define and find front time
locationAtFrontPassing = maximalNoise*NoiseFactor;
x_mns_x_tip = bsxfun(@minus, phediLocationMatrix, locationAtFrontPassing);
ToFindLeft = x_mns_x_tip;   ToFindLeft(x_mns_x_tip>0) = -inf;
ToFindRight = x_mns_x_tip;  ToFindRight(x_mns_x_tip<=0) = inf;
[~, closestFromLeft] = max(ToFindLeft,[],1);
[~, closestFromRight] = min(ToFindRight,[],1);
%--- use Thales law of triangles to interpolate to find time:
x_val_left = phediLocationMatrix(sub2ind(size(phediLocationMatrix), closestFromLeft, 1:length(closestFromLeft)));
x_val_right = phediLocationMatrix(sub2ind(size(phediLocationMatrix), closestFromRight, 1:length(closestFromRight)));
t_val_left = phediTimeVector(closestFromLeft);
t_val_right = phediTimeVector(closestFromRight);

tFrontPassing = mean(((t_val_right(:)-t_val_left(:))./(x_val_right(:)-x_val_left(:))).*(locationAtFrontPassing(:)-x_val_left(:))+t_val_left(:));


end