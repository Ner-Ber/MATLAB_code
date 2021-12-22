%---*** this script will plot the first and last rows of consecutive
%events. each event will be marked by a different color, first row by '*'
%and last by 's'

%% parameters
loadRows = 0;


%% load rows
if loadRows
    exper = my_dir;
    Exp = 5;
    EventsNums = 1:30;
    N = length(EventsNums);
    rowStorage = [];
    for i=1:N
        %     ROT = PhediStructCellUD{EventsNums(i)}.AsperityData.RowOverTime;
        CamMetaStruct = CameraMetaAllCams(exper{Exp});
        [preRowsTime, postRowsTime] = deal(5e-3);
        ROT = phantom_getRowOverTime(exper{Exp}, EventsNums(i), CamMetaStruct, preRowsTime, postRowsTime);
        first = ROT(1,:);
        last = ROT(end,:);
        newRows = [first(:)';last(:)'];
        rowStorage = cat(3,rowStorage,newRows);
    end
end

%% correlate among first rows
%--- correlate
allFirstRows = permute(rowStorage(1,:,:),[3 2 1]);
lagsVec = zeros(N,0);
dLagVec = zeros(N,0);
RowsCropped = cell(N,1);
spatialVecsCrop = cell(N,1);
for i = 2:N
    preRow = allFirstRows(i-1,:);
    postRow = allFirstRows(i,:);
    
    %--- use only high peaks:
    %-- create threshold
    Mpre = max(preRow);
    Mpost = max(postRow);
    higherLogicPre = preRow>Mpre/3;
    higherLogicPost = postRow>Mpost/3;
    %-- zero all below threshold
    preRowZ = preRow;   preRowZ(~higherLogicPre) = 0;
    postRowZ = postRow; postRowZ(~higherLogicPost) = 0;
    %-- get rid of sides below threshold:
    leftBoundPre = find(higherLogicPre,1,'first');
    rightBoundPre = find(higherLogicPre,1,'last');
    leftBoundPost = find(higherLogicPost,1,'first');
    rightBoundPost = find(higherLogicPost,1,'last');
    preRowCrp = preRowZ(leftBoundPre:rightBoundPre);
    postRowCrp = postRowZ(leftBoundPost:rightBoundPost);
    
    %--- create spatial vectors:
    XX = 1:length(preRow);
    XXpre = XX(leftBoundPre:rightBoundPre);
    XXpost = XX(leftBoundPost:rightBoundPost);
    
    %--- save the vectors
    RowsCropped{i-1} = preRowCrp;
    spatialVecsCrop{i-1} = XXpre;
    if i==N
        RowsCropped{i} = postRowCrp;
        spatialVecsCrop{i} = XXpost;
    end
    
    
    [r,lags] = xcorr(preRowCrp,postRowCrp);
    [~,I] = max(r);
    LAG = lags(I);
    lagsVec(i) = LAG;
    
    %--- consider initial difference
    dLag = XXpost(1)-XXpre(1);
    dLagVec(i) = dLag;
end

totalLag = lagsVec-dLagVec;
totalLagSUM = cumsum(totalLag);

%% plot shofted (hopefully collapsed)
relevantEvents = find(~cellfun(@isempty,PhediStructCellUD));
allPointsX = [];
allPointsY = [];
EvIdx = [];
phediSlope = [];
numPhedis = [];
numPhedisPlus = [];
numPhedisMinus = [];
EventColors = MyVaryColor(length(relevantEvents));
figure;
hold on;
for i=1:length(relevantEvents)
    %-- plot profile
    Xshifted = (1:L)+totalLagSUM(relevantEvents(i));
    plot(Xshifted,allFirstRows(relevantEvents(i),:),'Color',EventColors(i,:))
    
    %-- plot phedis
    currentIdxs = PhediStructCellUD{relevantEvents(i)}.PhediData.phedisInitialLocsInPix;
    currentSlope = PhediStructCellUD{relevantEvents(i)}.PhediData.slopeIncline;
    plot(Xshifted(currentIdxs),...
        allFirstRows(relevantEvents(i),currentIdxs),...
        'o','MarkerFaceColor',EventColors(i,:))
    
    EvIdx = cat(1,EvIdx,ones(length(currentIdxs),1)*relevantEvents(i));
    thisRoundSlope = PhediStructCellUD{relevantEvents(i)}.PhediData.slopeIncline';
    phediSlope = cat(1,phediSlope,thisRoundSlope);
    allPointsX = cat(1,allPointsX,Xshifted(currentIdxs)');
    allPointsY = cat(1,allPointsY,allFirstRows(relevantEvents(i),currentIdxs)');
    numPhedis = cat(1,numPhedis,length(currentIdxs));
    numPhedisPlus = cat(1,numPhedisPlus,nnz(thisRoundSlope==1));
    numPhedisMinus = cat(1,numPhedisMinus,nnz(thisRoundSlope==-1));
end



%% plot for clusters
% %-- cluster -1 slopes
% slopeLogical = phediSlope==-1;
% idx = kmeans(allPointsX(slopeLogical),max(numPhedisMinus));
% 
% totIdx = zeros(length(slopeLogical),1);
% totIdx(slopeLogical) = idx+100;
% 
% %-- cluster +1 slopes
% slopeLogical = phediSlope==1;
% idx = kmeans(allPointsX(slopeLogical),max(numPhedisPlus));
% totIdx(slopeLogical) = idx+200;
% 
% 
% relevantIdxs = unique(totIdx);
% relevantIdxs(relevantIdxs==0) = [];
% Cols = MyVaryColor(length(relevantIdxs),parula);
% for i=relevantIdxs(:)'
%     plot(allPointsX(totIdx==i),allPointsY(totIdx==i),'o','MarkerFaceColor',Cols(relevantIdxs==i,:));
% end
% 

%% cluster by find peaks
smthFac = 11;
smthSignalFirst = smooth(allFirstRows(1,:),smthFac);
smthSignalLast = smooth(allFirstRows(N,:),smthFac);
[pksMINfirst,locsMINfirst] = findpeaks(-smthSignalFirst);
[pksMAXfirst,locsMAXfirst] = findpeaks(smthSignalFirst);
[pksMINlast,locsMINlast] = findpeaks(-smthSignalLast);
[pksMAXlast,locsMAXlast] = findpeaks(smthSignalLast);

CloseMIN = interp1(locsMINlast,locsMINlast,locsMINfirst,'nearest');
meanValMIN = mean([CloseMIN(:),locsMINfirst(:)],2,'omitnan');

CloseMAX = interp1(locsMAXlast,locsMAXlast,locsMAXfirst,'nearest');
meanValMAX = mean([CloseMAX(:),locsMAXfirst(:)],2,'omitnan');

%--- create bins for sorting phedis
BinBorders = unique([0;meanValMAX(:);L+max(totalLagSUM)]);
Steps = multipleStepFunction(BinBorders);
XBins = Steps(allPointsX);


phedigroupsIndx = XBins.*phediSlope;

phedigroupsIndxUnq = unique(phedigroupsIndx);
Cols = MyVaryColor(length(phedigroupsIndxUnq),parula);
for i=phedigroupsIndxUnq(:)'
    plot(allPointsX(phedigroupsIndx==i),allPointsY(phedigroupsIndx==i),'o','MarkerFaceColor',Cols(phedigroupsIndxUnq==i,:));
end
