function UpdatedDataStruct = phedi_add_phediLocationVectors(DataStruct,varargin)
% UpdatedDataStruct = phedi_add_phediLocationVectors(DataStruct, CfSmooth)
%
% phedi_add_phediLocationVectors will ad x_mins_x_tip vectors for the phedi
% data of location and veclocity
%
% CfSmooth - number of pixel to smooth in (in big picture).

%% set defaults
CfSmooth = setDefaults4function(varargin,10);

%% load relevant data
%--- prep output
UpdatedDataStruct = DataStruct;
%--- load existing
BigPicRotStruct = DataStruct.BigPicRotStruct;
PhotoLocation = DataStruct.PhediData.PhotoLocation;
if isfield(DataStruct.PhediData,'t_mins_t_tip')
    t_mins_t_tip = DataStruct.PhediData.t_mins_t_tip;
else
    t_mins_t_tip = bsxfun(@minus,DataStruct.PhediData.timeVec(:),DataStruct.PhediData.t_tips(:)');
    UpdatedDataStruct = t_mins_t_tip;
end

%% time vector for each phedi
x_mins_x_tip = zeros(size(t_mins_t_tip));
relevantLogicals = zeros(size(t_mins_t_tip));
for i = 1:length(DataStruct.PhediData.t_tips)    
%     [x_min_Ttip_i, relevantLogical] = signal_ChangeTime2Space(t_mins_t_tip(:,i), PhotoLocation, BigPicRotStruct,CfSmooth);
    [x_min_Ttip_i, relevantLogical] = signal_ChangeTime2Space_mapping(t_mins_t_tip(:,i), PhotoLocation, BigPicRotStruct);
%     [x_min_Ttip_i, relevantLogical] = signal_ChangeTime2Space_mapping20181912(t_mins_t_tip(:,i), PhotoLocation, BigPicRotStruct);
    
    x_mins_x_tip(:,i) = x_min_Ttip_i(:);
    relevantLogicals(:,i) = relevantLogical(:);
end

%% add to structure
UpdatedDataStruct.PhediData.x_mins_x_tip = x_mins_x_tip;
UpdatedDataStruct.PhediData.inBlockLogicals = relevantLogicals;

%% update sgData;
UpdatedDataStruct = sg_updateSpaceVectors(UpdatedDataStruct,CfSmooth);
