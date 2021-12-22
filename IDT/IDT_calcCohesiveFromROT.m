function [cohesiveStruct] = IDT_calcCohesiveFromROT(RotStruct,region,varargin)
% [cohesiveStruct] = IDT_calcCohesiveFromROT(RotStruct,region,varargin)
%
% RotStruct - output of IDT_PlotRowOverTime.
% region - region of calculation needed, in meters
% IDT_calcCohesiveFromROT will calculate the cohesive zone
% VARARGIN = DEFAULT
% x_dA = -0.01  (the location where to spicify the exponent is at zero value)
% smt = 20. (smooth along frame. smt=1 for no smooth)


%% set defaults
[x_dA, smt] = setDefaults4function(varargin,-0.01, 20);

%% find x_tips
locationLogical = RotStruct.x<=max(region) & RotStruct.x>=min(region);
PixVector = 1:size(RotStruct.DataMat,2);
StepPixLogical = ismember(RotStruct.frontStepsPix,PixVector(locationLogical));
x_tipsPix = RotStruct.frontStepsPix(StepPixLogical);
x_tipsMeters = RotStruct.x(x_tipsPix);

% regionInPix = PixVector(locationLogical);
% StepLocLogical = RotStruct.frontStepsPix<=max(regionInPix) & RotStruct.frontStepsPix>=min(regionInPix);
% x_tipsPix = RotStruct.frontStepsPix(StepLocLogical);

%% get frames relevant for the region indicated

% locationLogical = RotStruct.x<=max(region) & RotStruct.x>=min(region);
% relevatFramesRegion = RotStruct.frontTime_interp(locationLogical);
% relevantFramesLogical = RotStruct.t<=(max(relevatFramesRegion)/RotStruct.fps)...
%     & RotStruct.t>=(min(relevatFramesRegion)/RotStruct.fps);
% relevantFrames = RotStruct.DataMatNorm(relevantFramesLogical,:);

RowNumVec = 1:size(RotStruct.DataMat,1);
RowNums = RotStruct.frontStepsTime(StepPixLogical);
relevantFramesLogical = ismember(RowNumVec,RowNums);
relevantFramesLogical = relevantFramesLogical(:);
relevantFrames = RotStruct.DataMatNorm(relevantFramesLogical,:);
%--- smooth frames
relevantFrames = conv2(relevantFrames,ones(1,smt)/smt,'same');

[N_frames,N_pixels] = size(relevantFrames);
relevantFramesLogical = relevantFramesLogical(:);

%% minus x tips
PixMatrix = repmat(PixVector,size(relevantFrames,1),1);
xVectorsPix = bsxfun(@minus,PixMatrix,x_tipsPix(:));
xVectorsMeters = xVectorsPix*RotStruct.res;

%% get values of delta_A
%--- define delta_A at:
[~,minIdx] = min(abs(xVectorsMeters-x_dA),[],2);
minLinIdx = sub2ind(size(xVectorsMeters),1:length(minIdx),minIdx(:)');
dA = 1-relevantFrames(minLinIdx);
dA = dA(:);
IdxOf_dARef = minLinIdx(:);

%% strech and find cohesive zone length
%--- strech the exponent to span from -1 to 0 by the following arithmetics
% ((relevantFrames-IdxOf_dARef)/(1-IdxOf_dARef))-1
renormalizedRows = bsxfun(@rdivide,...
                        bsxfun(@minus,...
                                relevantFrames,...
                                relevantFrames(IdxOf_dARef)...
                                ),...
                        (1-relevantFrames(IdxOf_dARef)))...
                                        -1;
%--- search for values where shifted normalized intensity equals -0.632%--- find relevant area to search in:
searchRegionLogic = xVectorsMeters<0 & xVectorsMeters>x_dA;
findXcRows = renormalizedRows-(-0.632);
relevantPixels = sum(searchRegionLogic(1,:));
[~,XcIndex] = min(abs(reshape(findXcRows(searchRegionLogic),[],relevantPixels)),[],2);
XcSingIndex = sub2ind([N_frames,relevantPixels],1:N_frames,XcIndex(:)');
xVectorsMetersRelevant = xVectorsMeters(searchRegionLogic);
Xc = abs(xVectorsMetersRelevant(XcSingIndex));

%% create structure for output
cohesiveStruct = struct;
[cohesiveStruct.Xc,...
    cohesiveStruct.dA,...
    cohesiveStruct.x_tipsMeters,...
    cohesiveStruct.relevantFramesLogical,...
    cohesiveStruct.xVectorsPix,...
    cohesiveStruct.xVectorsMeters,...
    cohesiveStruct.relevantFrames,...
    cohesiveStruct.renormalizedRows]...
 = deal(...
     Xc, dA, x_tipsMeters, relevantFramesLogical, xVectorsPix, xVectorsMeters, relevantFrames, renormalizedRows);

end