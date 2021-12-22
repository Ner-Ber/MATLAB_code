function [tFrontPassing , Cf_pixPerSec]= phedi_findFrontTimeFromCenterOfMass(...
    CM_TrajectoriesCrop,phediTimeVector,relevantFrames4Phedi,RowOverTime,phediInitialLocs)
%[tFrontPassing , Cf_pixPerSec]= phedi_findFrontTimeFromCenterOfMass(CM_Trajectories,phediTimeVector,relevantFrames4Phedi,RowOverTime,phediInitialLocs)
%
% phedi_findFrontTimeFromCenterOfMass will find the front passing time (in
% the time scale of the camera photographing the phedis) per phedi. 
% additionally, the function will output the speed (assuming constant) of
% the front when crossing the area photgraphed by the phedi-camera. 
%
% INPUTS:
% CM_Trajectories - trajectoris of the center of mass of the asperities.
%   this should be an mxn matrix, where m is the number of frames in
%   RowOverTime and n is the number of asperitis calculated. 
% phediTimeVector - time vector (in second) in the time scale of the camera
%   photographing the phedis. 
% relevantFrames4Phedi - the frames where the phedi movement was
%   calculated. 
% RowOverTime - 
% phediInitialLocs - initial locations of the phedis, in pix.
%
% OUTPUTS:
% tFrontPassing - a vector containing the time (in sec) at which the front
% passed the phedi. length(tFrontPassing)=number of phedis
% Cf_pixPerSec - front velocity in pix per sec.
%
%
% ASSUMPTIONS:
% 1) the speed and phedi crossing is calculated assuming that the drop in
% intensity in thr RowOverTime happens instantaniously as the front passes.
% 2) the front does not change its velocity over all the phedis. 
%

%% organize data
CM_TrajectoriesCrop = CM_TrajectoriesCrop(relevantFrames4Phedi,:);
CM_TrajectoriesNoNan = CM_TrajectoriesCrop(:,~isnan(sum(CM_TrajectoriesCrop)));
% CM_Traj_relevant = CM_TrajectoriesNoNan(:,[1,2,end-1,end]);
CM_Traj_relevant = CM_TrajectoriesNoNan;
%% create amplitude signals
CM_Traj_round = floor(CM_Traj_relevant);
pixCor = repmat(CM_Traj_round(1,:),size(CM_Traj_round,1),1);
% linCoor = sub2ind(size(RowOverTime),...
%     repmat(relevantFrames4Phedi(:),size(CM_Traj_round,2),1),pixCor(:));
linCoor = sub2ind(size(RowOverTime),...
    repmat(relevantFrames4Phedi(:),size(CM_Traj_round,2),1),CM_Traj_round(:));
AmpsSignal = reshape(RowOverTime(linCoor),length(relevantFrames4Phedi),[]);
AmpsSignal = bsxfun(@minus,AmpsSignal,mean(AmpsSignal(1:100,:)));
smoothFac = 11;
AmpsSignalSmooth = conv2(AmpsSignal,ones(smoothFac,1)/smoothFac,'same');

%% by falling below p2p value
refLength = 100;
fallingThresholds = 2*peak2peak(AmpsSignalSmooth(1:refLength,:),1);
fallingPoints = zeros(1,length(fallingThresholds));
for idx=1:length(fallingThresholds)
    y0 = -fallingThresholds(idx);
    y = AmpsSignalSmooth(:,idx);
    x = 1:length(AmpsSignalSmooth(:,idx));
    x0 = find_x0_at_y0(y0,y(:),x(:),'left');
    fallingPoints(idx) = x0;
end
%% find falling in signal strength by minimal derivative
AmpsSignalDeriv = conv2(AmpsSignalSmooth,[1 0 -1]','same');
border=20;
search_times = [-3 3]*1e-4;  % in [s]
searchLogical = phediTimeVector>=search_times(1) & phediTimeVector<=search_times(2);
% [~, minAmpsSignalDeriv] = min(AmpsSignalDeriv(border:end-border,:),[],1);
AmpsSignalDeriv2search = AmpsSignalDeriv;
AmpsSignalDeriv2search(~searchLogical,:) = inf;
AmpsSignalDeriv2search([1:border, (end-border):end],:) = inf;
[~, minAmpsSignalDeriv] = min(AmpsSignalDeriv2search,[],1);
% minAmpsSignalDeriv= minAmpsSignalDeriv+border;



%% assume constant velocity and assign time for each asperity
boundary1 = mean(CM_Traj_relevant(1,1:2));
boundary1_frame = round(mean(minAmpsSignalDeriv(1:2)));
boundary1_time = interp1(1:length(phediTimeVector),phediTimeVector,boundary1_frame);
boundary2 = mean(CM_Traj_relevant(1,3:4));
boundary2_frame = round(mean(minAmpsSignalDeriv(3:4)));
boundary2_time = interp1(1:length(phediTimeVector),phediTimeVector,boundary2_frame);
%--- create front time passing for each pixel:
tFrontPassing = interp1([boundary1 boundary2], [boundary1_time,boundary2_time], phediInitialLocs,'linear','extrap');
Cf_pixPerSec = (phediInitialLocs(end)-phediInitialLocs(1))/(tFrontPassing(end)-tFrontPassing(1));


end