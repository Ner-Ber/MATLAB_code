function StatStruct = timeSamplesGetStatistics(timeSamples, samples2Use)
% timeSamples - is a NxM matrix contaning M time series measurments, each
%       measurmant in the length of N time steps. 
% samples2Use - is a vector of indicies implying what samples (columns from
%       'timeSamples') to use.


if nargin<2
    samples2Use= 1:size(timeSamples,2);
end
RelvntTimeSam = timeSamples(:,samples2Use);

%--- calc mean signal
meanSig = mean(RelvntTimeSam,2);
StatStruct.mean = meanSig;

%-- calc std of measurments
stdSig = std(RelvntTimeSam,0,2);
StatStruct.std= stdSig;

%-- calc minimal trajectory of measurments
minSig = min(RelvntTimeSam,[],2);
StatStruct.min= minSig;

%-- calc maximal trajectory of measurments
maxSig = max(RelvntTimeSam,[],2);
StatStruct.max = maxSig;





