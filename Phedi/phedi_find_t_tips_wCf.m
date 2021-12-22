function [intersectionFrame, crossTimeInFrm, crossLocInPix] = phedi_find_t_tips_wCf(RowOverTime,CfPixPerFrame,phediInitialLocs,PhediLocationPix,phediFrames)
% [intersectionFrame, crossTimeInFrm, crossLocInPix] = phedi_find_t_tips_wCf(RowOverTime,CfPixPerFrame,phediInitialLocs,PhediLocationPix,phediFrames)
%
% phedi_find_t_tips_wCf will find the intersection frame of the front in
% the small RowOverTime and find t_tips and their location in units of
% frames.


[CfPixPerFrm, frontAppearenceFrame, X, Y] = phantom_findCf(RowOverTime,1:size(RowOverTime,1));
%--- model of whoch to find the frame enetering:
model = @(InterFrame,x) bsxfun(@plus,InterFrame,(1./CfPixPerFrame).*x);
%--- define entering frames to check:
Dy = ceil(size(RowOverTime,2)./CfPixPerFrame);
% InterFrameVec = linspace(-Dy,(size(RowOverTime,1)+Dy),1e4);
InterFrameVec = -Dy:0.1:(size(RowOverTime,1)+Dy);
%--- find time with highest density of dots suspected as front
[N,edges] = histcounts(Y,round(size(RowOverTime,1)/100));
[~,maxNidx] = max(N);
borders = [edges(maxNidx),edges(maxNidx+1)];
border_expand = [borders(1)-diff(borders), borders(2)+diff(borders)];
inRange = Y>min(border_expand) & Y<max(border_expand);
Yin = Y(inRange);
Xin = X(inRange);
% P0 = mean(InterFrameVec);
% [lb,ub] = deal(min(InterFrameVec),max(InterFrameVec));
% model = @(P,x) P(1)+(1./CfPixPerFrame).*x;
% opts = optimset('Display','off');
% P = lsqcurvefit(model,P0,X,Y,lb,ub,opts);


% modelsDataDiff = bsxfun(@minus,Y(:)',model(InterFrameVec(:),X(:)'));
% STD = std(abs(modelsDataDiff),0,2);
% % logicalSTD = modelsDataDiff>repmat(STD,1,size(modelsDataDiff,2));
% logicalSTD = abs(modelsDataDiff)>100;
% modelsDataDiff_selected = modelsDataDiff;
% modelsDataDiff_selected(logicalSTD) = nan;
% %--- rms ignoring nans
% numNans = sum(isnan(modelsDataDiff_selected),2);
% N_differential = size(logicalSTD,2)-numNans;
% S = sum(modelsDataDiff_selected.^2,2,'omitnan');
% R = sqrt(S./N_differential);
% RMS = rms(modelsDataDiff_selected,2);
% %--- find minimal RMS as indicator for best fit:
% [~,I] = min(R);
modelsDataDiff = bsxfun(@minus,Yin(:)',model(InterFrameVec(:),Xin(:)'));
RMS = rms(modelsDataDiff,2);
%--- find minimal RMS as indicator for best fit:
[~,I] = min(RMS);
intersectionFrame = round(InterFrameVec(I));

Front_t = @(x) intersectionFrame+(1./CfPixPerFrame).*x;

%% find t-t_tips with new model
PhediLocPixOnROT = bsxfun(@plus, PhediLocationPix, phediInitialLocs(:)');
phediTimes = repmat(phediFrames(:),1,length(phediInitialLocs));
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


end