function NewDataStruct = phedi_addSmoothedFunctions2Struct(DataStruct,varargin)
%phedi_addSmoothedFunctions2Struct will add to DataStruct (the output of 
% Movie_phedi_from_folder_2_data) the phedi location and velocity after
% smoothenong phediLocation by a known factor. 

%% set defaults
[Nsmt] = setDefaults4function(varargin,4);

%% backup original data
NewDataStruct = DataStruct;
NewDataStruct.PhediDataPreSmth = DataStruct.PhediData;

%% get data from current structure
PhediLocation = NewDataStruct.PhediData.PhediLocation;
PhediLocationPix = NewDataStruct.PhediData.PhediLocationPix;

%% smooth phedi location 
PhediLocation = conv2(PhediLocation,ones(Nsmt,1)/Nsmt,'same');
PhediLocation([1:Nsmt,end-Nsmt:end],:) = nan;
NewDataStruct.PhediData.PhediLocation = PhediLocation;

PhediLocationPix = conv2(PhediLocationPix,ones(Nsmt,1)/Nsmt,'same');
PhediLocationPix([1:Nsmt,end-Nsmt:end],:) = nan;
NewDataStruct.PhediData.PhediLocationPix = PhediLocationPix;

%% add new phedi velocity
% [PhediVelocity, timeVec_4vel] = phedi_calcVelocity(PhediLocation,timeVec);
NewDataStruct = phedi_add_velocity_to_struct(NewDataStruct);

end