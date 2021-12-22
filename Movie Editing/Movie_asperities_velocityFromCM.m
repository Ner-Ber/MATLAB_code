function [V_cm, t_mins_t_tip, x_min_x_tip] = Movie_asperities_velocityFromCM(DataStruct,varargin)
    
    %% set defaults
    [NsmthLoc] = setDefaults4function(varargin,3);
    
    %% load relevant data
    center_of_mass = DataStruct.AsperityData.center_of_mass;
    res = DataStruct.ExperimentData.res;
    CM_time = DataStruct.AsperityData.timeCount;
    spatialVec = DataStruct.AsperityData.spatialVec;
    timeCount = DataStruct.AsperityData.timeCount;
    center_of_mass_loc = center_of_mass./res;
    BigPicRotStruct = DataStruct.BigPicRotStruct;
    PhotoLocation = DataStruct.PhediData.PhotoLocation;
    Ncm = size(center_of_mass_loc,2);
    %     PhediLocation = DataStruct.PhediData.PhediLocation;
    %     measuredPhedisFromPlot = DataStruct.PhediData.measuredPhedisFromPlot;
    %% calc velocity
    %-- -smooth location
    center_of_mass_locSMTH = conv2(center_of_mass_loc,ones(NsmthLoc,1)./NsmthLoc,'same');
    
    
    CM_timeMAT = repmat(CM_time,1,Ncm);
    V_cm = diff(center_of_mass_locSMTH)./diff(CM_timeMAT);
    CM_time_4vel = movmean(CM_time,2,'Endpoints','discard');
    
    %% find t_tips
    [x_tips,t_tips] = phedi_find_t_tips_wPhotoLocation(BigPicRotStruct,PhotoLocation,spatialVec,...
        bsxfun(@minus,center_of_mass_locSMTH,center_of_mass_locSMTH(1,:)),center_of_mass_locSMTH(1,:),...
        timeCount);
    %% create time vecs for each CM
    CM_time_4velMAT = repmat(CM_time_4vel,1,size(center_of_mass_locSMTH,2));
    t_mins_t_tip = bsxfun(@minus,CM_time_4velMAT,t_tips(:)');
    
    %% x_min_x_tip
    x_min_x_tip = [];
    for i=1:Ncm
        x_min_Ttip = signal_ChangeTime2Space_mapping(t_mins_t_tip(:,i), PhotoLocation, BigPicRotStruct);
        x_min_x_tip = cat(2,x_min_x_tip,x_min_Ttip);
    end
end