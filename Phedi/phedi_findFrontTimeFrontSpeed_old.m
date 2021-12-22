function [tFrontPassing , CfPixPerFrm, frontAppearenceFrame]= phedi_findFrontTimeFrontSpeed(RowOverTime,relevantFrames4Phedi,phediInitialLocs,PhediLocationPix)
% [tFrontPassing , Cf_pixPerFrm]= phedi_findFrontTimeFrontSpeed(RowOverTime,relevantFrames4Phedi,phediInitialLocs)
%
% !!!! all inputs are in units of pixels and frames !!!!


[CfPixPerFrm, frontAppearenceFrame, X, Y] = phantom_findCf(RowOverTime,relevantFrames4Phedi);
%% equation of front time:
Front_t = @(x) frontAppearenceFrame+(1/CfPixPerFrm).*x;
%% find where phedi trajectories cross the equation of the front
PhediLocPixOnROT = bsxfun(@plus, PhediLocationPix, phediInitialLocs(:)');
phediTimes = repmat(relevantFrames4Phedi(:),1,length(phediInitialLocs));
equationTimes = Front_t(PhediLocPixOnROT);
%--- closest measurement of phedi and front
diffTime = equationTimes-phediTimes;
%--- prepare a vector to find from it crossings:
diffTimeLogicDoub = double(diffTime>0)*2-1;
diffTimeStep = conv2(diffTimeLogicDoub,[-1 1]'/2,'same');
diffTimeStep([1:2,end-1:end],:)=0;
crossFrameCell = cell(1,length(phediInitialLocs));
%--- find frame where phedi trajectories are crossing front:
for i=1:length(phediInitialLocs); crossFrameCell{i} = find(diffTimeStep(:,i)); end;
%-- assume that if there is more than one crossing, the trajectory is
%negligabily wiggily and the crossings are close:
crossFramePartial = round(cellfun(@mean,crossFrameCell));
crossFrame = [crossFramePartial;crossFramePartial+1];
columnIdx = [1:length(phediInitialLocs);1:length(phediInitialLocs)];
LinIdx = sub2ind(size(PhediLocPixOnROT),crossFrame(:),columnIdx(:));
if nnz(isnan(LinIdx))>0
    if nnz(isnan(LinIdx))>round(length(LinIdx)/1.5)
        [tFrontPassing , CfPixPerFrm, frontAppearenceFrame] = deal(0);
        return;
    end
    LinIdx(isnan(LinIdx)) = 1;
end
crossLoc = PhediLocPixOnROT(LinIdx);
crossLocReshp = reshape(crossLoc,2,[]);
crossFramesRot = relevantFrames4Phedi(crossFrame);
%% create linear equations for crossing points
slopes = (crossFramesRot(2,:)-crossFramesRot(1,:))./(crossLocReshp(2,:)-crossLocReshp(1,:));
n_s = -slopes.*crossLocReshp(1,:)+crossFramesRot(1,:);
%% find crossing points of linear interpolation and front location:
locs_0 = (n_s-frontAppearenceFrame)./((1/CfPixPerFrm)-slopes);
frames_0 = slopes.*locs_0+n_s;
%% assign to crossing points
tFrontPassing = frames_0;
%% eliminate nans from infinite slope
nanLocs = isnan(frames_0);
tFrontPassing(nanLocs) = Front_t(crossLocReshp(1,nanLocs));


% 
% %% get interpolated crossing points
% %--- linear line interpolation of the phedi trajectory crossing:
% tFrontPassing = Front_t(phediInitialLocs);
