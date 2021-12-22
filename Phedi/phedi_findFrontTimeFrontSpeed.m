function [tFrontPassing , crossLocInPix, CfPixPerFrm, frontAppearenceFrame] = phedi_findFrontTimeFrontSpeed(RowOverTime,relevantFrames4Phedi,phediInitialLocs,PhediLocationPix)
% [tFrontPassing , CfPixPerFrm, frontAppearenceFrame]= phedi_findFrontTimeFrontSpeed(RowOverTime,relevantFrames4Phedi,phediInitialLocs,PhediLocationPix)
%
% !!!! all inputs are in units of pixels and frames !!!!


% [CfPixPerFrm, frontAppearenceFrame, X, Y] = phantom_findCf(RowOverTime,relevantFrames4Phedi);
[CfPixPerFrm, frontAppearenceFrame, X, Y] = phantom_findCf(RowOverTime,1:size(RowOverTime,1));
%% equation of front time:
Front_t = @(x) frontAppearenceFrame+(1/CfPixPerFrm).*x;
%% find where phedi trajectories cross the equation of the front
PhediLocPixOnROT = bsxfun(@plus, PhediLocationPix, phediInitialLocs(:)');
phediTimes = repmat(relevantFrames4Phedi(:),1,length(phediInitialLocs));
equationTimes = Front_t(PhediLocPixOnROT);
%--- find intersections of front and phedi trajectory:
space_vec = linspace(-1,(size(RowOverTime,2)+1),25);
crossTimeInFrm = [];
crossLocInPix = [];
for i=1:length(phediInitialLocs)
    [x1, t1, x2, t2] = deal(space_vec,Front_t(space_vec),PhediLocPixOnROT(:,i),phediTimes(:,i));
    [x0,t0] = intersections(x1,t1,x2,t2);
    %-- assume that if there is more than one crossing, the trajectory is
    %negligabily wiggily and the crossings are close:
    crossTimeInFrm = cat(1,crossTimeInFrm,mean(t0));
    crossLocInPix = cat(1,crossLocInPix,mean(x0));
end
%% deal with NaNs
%--- nans will appear when the phedi trajectory and the front don't
%intersect. in this case they will be replace by a linear exptrapolation
%just in order to put a number there. THIS IS VERY NOT ACCURATE. 

nanLocations = isnan(crossTimeInFrm);
%--- compare front time at each phedi to the mean time of phedis
meanPhediTime = mean(relevantFrames4Phedi);
differenceSign = sign(meanPhediTime-Front_t(phediInitialLocs));    % when negative front is below phedi
tFrontPassing = crossTimeInFrm;
FrontBelow = nanLocations & differenceSign==-1;
FrontAbove = nanLocations & differenceSign==1;
tFrontPassing(FrontBelow) = Front_t(PhediLocPixOnROT(end,FrontBelow));
tFrontPassing(FrontAbove) = Front_t(PhediLocPixOnROT(1,FrontAbove));

