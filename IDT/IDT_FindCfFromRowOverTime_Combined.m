function frontStructure = IDT_FindCfFromRowOverTime(RowOverTime,smt)
%frontStructure = IDT_FindCfFromRowOverTime(RowOverTime,smt)
%IDT_FindCfFromRowOverTime - will return a structure containing data about
%the front:
% fronLoc,StepsPix,StepTime,FrontVel,VelLoc.
% all units are in pixles and frames

if (nargin<2)                                    %% default is smooth factor 4
    smt=4;
end
S = size(RowOverTime);

%% search for burnt frames and replace with nans
%-- check on a few columns (pixels) if there is any need at all:
[~,~,~,p] = findpeaks(RowOverTime(:,round(S(2)/2)));
%--- if it is suspected, check and replace row with nans
if max(p)/mean(p)>=10
    [~,~,~,pMat] = MyFindpeaks(RowOverTime(:,round(S(2)/3):round(2*S(2)/3)));
    [~, I_max_p] = max(pMat,[],1,'omitnan');
    [M,F] = mode(I_max_p);
    if F/length(round(S(2)/3):round(2*S(2)/3))>0.75
        RowOverTime((M-1):(M+1),:) = nan;
    end
end



%% Find approximate front area by using first derivative

%--- use normalized ROT to find first derivative
RowOverTimeNorm=bsxfun(@rdivide,RowOverTime,RowOverTime(1,:));
TipThresh = 1.2;
RowOverTimeNorm(RowOverTimeNorm>TipThresh) = TipThresh;
RowOverTimeNorm(RowOverTimeNorm<0) = 0;

%--- first derivative of normalized
FrstDrvROT=conv2(RowOverTimeNorm,[1 0 -1]','same');

%--smooth the picture
SmoothedROT=conv2(RowOverTime,ones(smt)/smt^2,'same');

%--2nd derivative
ScndDrvROT=conv2(SmoothedROT,[1 -2 1]','same');

%--- eliminate
n = 5; %remove from edge (since there are some edge effects from convolution)
ScndDrvROT(:,1:n) = nan; ScndDrvROT(:,(end-n+1):end) = nan; ScndDrvROT(1:n,:) = nan; ScndDrvROT((end-n+1):end,:) = nan;
FrstDrvROT(:,1:n) = nan; FrstDrvROT(:,(end-n+1):end) = nan; FrstDrvROT(1:n,:) = nan; FrstDrvROT((end-n+1):end,:) = nan;

%--- iterate until finding real front area usinf forst derivative (not frame slip)
iter = 1;
realFront=0;
while realFront==0 || iter<100
    [~,ImnsFirst]=min(FrstDrvROT,[],1,'omitnan');
    fronLocByMin = ImnsFirst;
    
    %--- is found front actual front or frame slip?
    %-- check this on central third of row
    Nthird = round(size(RowOverTime,2)/3);
    fronDiff = diff(fronLocByMin);
    fronDiffThird = fronDiff(Nthird:2*Nthird);
    
    if nnz(fronDiffThird)/length(fronDiffThird)>0.01
        realFront=1;
        break
    end
    FrstDrvROT(ImnsFirst,1:size(FrstDrvROT,2)) = nan;
    %--- promote counter
    iter = iter+1;
end

%% use approximate area of front to search in
%--- crop the area of the front fount by first derivative
FindFrontInterval = 10;
IntervalVec = (-FindFrontInterval:FindFrontInterval)';
ShiftsMat = repmat(IntervalVec,1,size(FrstDrvROT,2));
% fronLocByMin = smooth(fronLocByMin,5);
fronLocByMinMat = repmat(fronLocByMin(:)',2*FindFrontInterval+1,1);
fronLocByMinMat = fronLocByMinMat-ShiftsMat;   % row coordinates
columnCoordinateVec = repmat(1:size(FrstDrvROT,2),2*FindFrontInterval+1,1); % column coordinates
%--- if indecies are larger than ROT:
fronLocByMinMat(fronLocByMinMat<1) = 1;
fronLocByMinMat(fronLocByMinMat>size(FrstDrvROT,1)) = size(FrstDrvROT,1);
columnCoordinateVec(columnCoordinateVec<1) = 1;
columnCoordinateVec(columnCoordinateVec>size(FrstDrvROT,2)) = size(FrstDrvROT,1);

linearInd = sub2ind(size(ScndDrvROT), fronLocByMinMat(:), columnCoordinateVec(:));
ScndDrvROTCropped = ScndDrvROT(linearInd);
ScndDrvROTCropped = reshape(ScndDrvROTCropped,2*FindFrontInterval+1,size(FrstDrvROT,2));

%% find exact location of front using second derivative (zero crossing eadge detection)
%--- find fron location by using second derivative
[~,Imxs]=max(ScndDrvROTCropped,[],1,'omitnan');
[~,Imns]=min(ScndDrvROTCropped,[],1,'omitnan');
%--- find relevant location considering shifts
ImxsFixed = fronLocByMinMat(sub2ind(size(ShiftsMat), Imxs, 1:size(ShiftsMat,2)));
ImnsFixed = fronLocByMinMat(sub2ind(size(ShiftsMat), Imns, 1:size(ShiftsMat,2)));

%--find midpoint
fronLoc=floor((ImxsFixed+ImnsFixed)./2);    %this produces a "stairs" graph


%--take only the startpoint of the stair
StepsPix=find(diff(fronLoc));
StepTime=fronLoc(StepsPix)+1;   % the +1 is so the frame will coincide with the image

%% compuste front velocity
%--compute velocity in pixels/frame
dX=diff(StepsPix);
dT=diff(StepTime);
FrontVel=dX./dT;
VelLoc=movmean(StepsPix,2,'Endpoints','discard');

frontStructure = struct;
[frontStructure.fronLoc,frontStructure.StepsPix,frontStructure.StepTime,frontStructure.FrontVel,frontStructure.VelLoc] =...
    deal(fronLoc,StepsPix,StepTime,FrontVel,VelLoc);
end