function [cohesiveStruct] = IDT_calcCohesiveFromROT_time2space(RotStruct,region,varargin)
% [cohesiveStruct] = IDT_calcCohesiveFromROT_time2space(RotStruct,region,varargin)
%
% RotStruct - output of IDT_PlotRowOverTime.
% region - region of calculation needed, in meters
% IDT_calcCohesiveFromROT will calculate the cohesive zone
% VARARGIN = DEFAULT
% x_dA = -0.01  (the location where to spicify the exponent is at zero value)
% smt = 20. (smooth along frame. smt=1 for no smooth)


%% set defaults
[x_dA, smt] = setDefaults4function(varargin,-0.01, 20);

%% get relevant frames
locationLogical = RotStruct.x<=max(region) & RotStruct.x>=min(region);
relevantFrames = RotStruct.DataMatNorm(:,locationLogical)';
%--- smooth frames
relevantFrames = conv2(relevantFrames,ones(1,smt)/smt,'same');

%% get matching time vectors
t_tipsSeconds = RotStruct.frontTime_interp(locationLogical)/RotStruct.fps;
timeMatrix = repmat(RotStruct.t(:)',size(relevantFrames,1),1);
tVectors = bsxfun(@minus,timeMatrix,t_tipsSeconds(:));

%% get matching space vectors
%--- interpolate Cf vector for corrosponding space
Cf = RotStruct.frontVel_interp*RotStruct.res*RotStruct.fps;
Cf_interp = interp1(1.5:1:(length(Cf)+0.5),Cf,1:length(locationLogical));
relevantCf = Cf_interp(locationLogical);
xVectorsMeters = bsxfun(@times, tVectors,relevantCf(:));

x_tipsMeters = t_tipsSeconds.*relevantCf;
[N_frames,N_pixels] = size(relevantFrames);

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
    cohesiveStruct.xVectorsMeters,...
    cohesiveStruct.relevantFrames,...
    cohesiveStruct.renormalizedRows]...
 = deal(...
     Xc, dA, x_tipsMeters, xVectorsMeters, relevantFrames, renormalizedRows);

end