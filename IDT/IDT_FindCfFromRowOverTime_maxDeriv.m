function frontStructure = IDT_FindCfFromRowOverTime(RowOverTime,smt)
%frontStructure = IDT_FindCfFromRowOverTime(RowOverTime,smt)
%
% Method: in each row (frame) in the normalized RowOverTime find the points
% where the value changes from 1 to something lower. 
% specifically, since there is noise, the row is smoothened (Nsmooth
% variable in the code determines how smooth it will be) and then rounded
% up to 1 decimal. effectively the code will find the transition between
% 0.95 to somthing lower.
%
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

%% find front by finding falling from values of 1
%--- get ranges to find the step in. take each point on the front and
%search for the change in value around it:
range2search = 8;
range2searchMAT = -range2search:range2search;
rangesMAT = [1:size(RowOverTime,2)]';
ColumnsIdx = bsxfun(@plus,repmat(rangesMAT,1,length(range2searchMAT)),range2searchMAT);
ColumnsIdx(ColumnsIdx<1) = 1;
ColumnsIdx(ColumnsIdx>size(RowOverTime,2)) = size(RowOverTime,2);
RowsIdx = repmat(ImnsFirst(:),1,length(range2searchMAT));
LinInd = sub2ind(size(RowOverTime),RowsIdx(:),ColumnsIdx(:));

%--- create a ROT compatible for this type of search
Nsmooth = 11;
RowOverTimeNormSmooth = conv2(RowOverTimeNorm,ones(1,Nsmooth)/Nsmooth,'same');
RowOverTimeNormRound = round(RowOverTimeNormSmooth,1);
RowOverTimeNormRound(isnan(FrstDrvROT)) = nan;

%--- find steps by deriving the rounded ROT
RowOverTimeNormRoundLogic = RowOverTimeNormRound;
RowOverTimeNormRoundLogic(RowOverTimeNormRoundLogic~=1) = 0;
mask = [-1 1];
ROT_LogicDeriv = conv2(RowOverTimeNormRoundLogic,mask,'same');
[maxVal, maxIdx] = max(abs(ROT_LogicDeriv),[],1);
StepLogical = maxVal~=0;
xx = 1:S(2);
StepRow = maxIdx(StepLogical);
StepColTMP = xx(StepLogical);
%--- fix step location by finding closest 1 value in row
StepColFILL = zeros(size(StepColTMP));
for i=1:length(StepColTMP)
    logicalRow = RowOverTimeNormRoundLogic(StepRow(i),:);
    relevantDistances = xx;
    relevantDistances(~logicalRow)=inf;
    [~,StepColFILL(i)] = min(abs(relevantDistances-StepColTMP(i)));
end

%--- find actuatl pixel in frames
[StepCol,ia,~] = unique(StepColFILL);
StepsPix = StepCol;
StepTime = StepRow(ia);


%% compuste front velocity
%--compute velocity in pixels/frame
dX=diff(StepsPix);
dT=diff(StepTime);
FrontVel=dX./dT;
VelLoc=movmean(StepsPix,2,'Endpoints','discard');

frontStructure = struct;
[frontStructure.fronLoc,frontStructure.StepsPix,frontStructure.StepTime,frontStructure.FrontVel,frontStructure.VelLoc] =...
    deal(ImnsFirst,StepsPix,StepTime,FrontVel,VelLoc);
end