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

%% edit dtuff on row over time
%--smooth the picture
SmoothedROT=conv2(RowOverTime,ones(smt)/smt^2,'same');

%--2nd derivative
ScndDrvROT=conv2(SmoothedROT,[1 -2 1]','same');

n = 5; %remove from edge (since there are some edge effects from convolution)
ScndDrvROT(:,1:n) = nan; ScndDrvROT(:,(end-n+1):end) = nan; ScndDrvROT(1:n,:) = nan; ScndDrvROT((end-n+1):end,:) = nan;

%--- iterate until finding real fron (not frame slip)
iter = 1;
realFront=0;
while realFront==0 || iter>100
    [~,Imxs]=max(ScndDrvROT,[],1,'omitnan');
    [~,Imns]=min(ScndDrvROT,[],1,'omitnan');
    
    %--find midpoint
    fronLoc=floor((Imxs+Imns)./2);    %this produces a "stairs" graph
    
    %--- is found front actual front or frame slip?
    %-- check this on central third of row
    Nthird = round(size(RowOverTime,2)/3);
    fronDiff = diff(fronLoc);
    fronDiffThird = fronDiff(Nthird:2*Nthird);
    
    %--- debbuging
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
    ScndDrvROT(makeNaNs,[1:size(ScndDrvROT,2),1:size(ScndDrvROT,2)]) = nan;
    %--- promote counter
    iter = iter+1;
end
%--take only the startpoint of the stair
StepsPix=find(diff(fronLoc));
StepTime=fronLoc(StepsPix)+1;   % the +1 is so the frame will coincide with the image

%--compute velocity in pixels/frame
dX=diff(StepsPix);
dT=diff(StepTime);
FrontVel=dX./dT;
VelLoc=movmean(StepsPix,2,'Endpoints','discard');

frontStructure = struct;
[frontStructure.fronLoc,frontStructure.StepsPix,frontStructure.StepTime,frontStructure.FrontVel,frontStructure.VelLoc] =...
    deal(fronLoc,StepsPix,StepTime,FrontVel,VelLoc);
end