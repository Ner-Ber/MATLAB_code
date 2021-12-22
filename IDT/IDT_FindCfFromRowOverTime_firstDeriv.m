function frontStructure = IDT_FindCfFromRowOverTime(RowOverTime,smt)
%frontStructure = IDT_FindCfFromRowOverTime(RowOverTime,smt)
%IDT_FindCfFromRowOverTime - will return a structure containing data about
%the front:
% fronLoc,StepsPix,StepTime,FrontVel,VelLoc.
% all units are in pixles and frames

if (nargin<2)                                    %% default is smooth factor 4
    smt=4;
end

% RowOverTime=conv2(RowOverTime,ones(smt)/smt^2,'same');

%--- use normalized ROT
RowOverTime=bsxfun(@rdivide,RowOverTime,RowOverTime(1,:));
TipThresh = 1.2;
RowOverTime(RowOverTime>TipThresh) = TipThresh;
RowOverTime(RowOverTime<0) = 0;

%--smooth the picture
SmoothedROT=conv2(RowOverTime,ones(smt)/smt^2,'same');
% SmoothedROT = RowOverTime;

%--2nd derivative
FrstDrvROT=conv2(SmoothedROT,[1 0 -1]','same');

FrstDrvROT(:,1:3) = nan; FrstDrvROT(:,(end-2):end) = nan; FrstDrvROT(1:3,:) = nan; FrstDrvROT((end-2):end,:) = nan;

n = 5; %remove from edge (since there are some edge effects from convolution)
%--- iterate until finding real fron (not frame slip)
iter = 1;
realFront=0;
while realFront==0 || iter>100
%     [~,Imxs]=max(FrstDrvROT(n:end-n,:),[],1,'omitnan');
    [~,Imns]=min(FrstDrvROT(n:end-n,:),[],1,'omitnan');
%     Imxs = Imxs+n-1;
    Imns = Imns+n-1;
    
    %--find midpoint
    fronLoc=floor(Imns);    %this produces a "stairs" graph
    
    %--- is found front actual front or frame slip?
    %-- check this on central third of row
    Nthird = round(size(RowOverTime,2)/3);
    fronDiff = diff(fronLoc);
    fronDiffThird = fronDiff(Nthird:2*Nthird);
    
%     %--- debbuging
%     figure; imagesc(ScndDrvROT);
%     hold on;
%     plot(1:1920,Imxs,'r*');
%     plot(1:1920,Imns,'g*');
%     disp('real front criterea'); nnz(fronDiffThird)/length(fronDiffThird)>0.01
    
    
    if nnz(fronDiffThird)/length(fronDiffThird)>0.01
        realFront=1;
        break
    end
    makeNaNs = [Imxs,Imns];
    FrstDrvROT(makeNaNs,[1:size(FrstDrvROT,2),1:size(FrstDrvROT,2)]) = nan;
    %--- promote counter
    iter = iter+1;
end
%--take only the startpoint of the stair
StepsPix=find(diff(fronLoc));
StepTime=fronLoc(StepsPix)+1;   % the +1 is so the frame will coincide with the image
% StepTime=fronLoc(StepsPix);

%--compute velocity in pixels/frame
dX=diff(StepsPix);
dT=diff(StepTime);
FrontVel=dX./dT;
VelLoc=movmean(StepsPix,2,'Endpoints','discard');

frontStructure = struct;
[frontStructure.fronLoc,frontStructure.StepsPix,frontStructure.StepTime,frontStructure.FrontVel,frontStructure.VelLoc] =...
    deal(fronLoc,StepsPix,StepTime,FrontVel,VelLoc);
end