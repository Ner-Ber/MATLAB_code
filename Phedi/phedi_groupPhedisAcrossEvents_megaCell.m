function [phediGroupingCell, slopesByOrder, GroupsOrder, phediIntlLocByEv] = phedi_groupPhedisAcrossEvents_megaCell(PhediStructCellRelevant, varargin)
    %phediGroupingCell = phedi_groupPhedisAcrossEvents(experNum,PhediStructCellOfCell, varargin)
    %
    %phedi_groupPhedisAcrossEvents will group identicle phedis across different
    %events in the same experiment.
    %
    % THIS FUNCTION WORLS ONLY IF THE MULTIPLE CELL ARRAYS WERE PHOTOGRAPHED
    % WITH THE SAME PREFERENCES (fps, resolution...)
    %
    % phedi_groupPhedisAcrossEvents reads RowOverTime of the enlargment optical
    % system from scratch inorder to avoid smoothing if was done. This is
    % needed for a corelation check between images of different events.
    %
    % PhediStructCellOfCell - is a cell array containing multiple cell arrays
    % which are the known DataStruct (resulting from the function
    % Movie_phedi_from_folder_2_data_parameterStruct) where each structure of
    % an event is located in the corrosponding cell number.
    %
    % varargin - contains an (optional):
    % EventsNums - event vector to consider. default is 1:30.
    % rowsFlag - load 'first' rows, 'last' rows or 'all'.
    %
    % *** this function may acidentally separate phedis to different groups,
    % but (most likely) will not mix them.
    % additionally, it may ignore some phedis as a matter of caution (if it
    % recognized more than a single phedi of the same group in the same
    % event***
    %
    %
    %
    %   THIS IS AN UPDATED VERSION (2018-05-19):
    %   INSTEAD OF TAKING A 'CELL OF CELLS' IT WILL TAKE A SORTED CELL (A CELL
    %   CONTAINING ONLY RELEVANT EVENTS).
    %
    
    
    %% defaults
    [EventsNums, rowsFlag] = setDefaults4function(varargin,1:30,'first');
    %% load rows
    N = length(PhediStructCellRelevant);
    rowStorage = [];
    
    %-- get camera data
    exper_dirCell = {};
    eventNumCell = {};
    try
        CamMetaStruct = PhediStructCellRelevant{1}.CamMeta;
    catch
        CamMetaStruct = CameraMetaAllCams(PhediStructCellRelevant{1}.ExperimentData.ExpHour);
    end
    
    for e = 1:N
        exper_dirCell_e = PhediStructCellRelevant{e}.ExperimentData.ExpHour;
        exper_dirCell{e} = exper_dirCell_e;
        eventNumCell_e = PhediStructCellRelevant{e}.ExperimentData.eventNum;
        eventNumCell{e} = eventNumCell_e;
        try
            if strcmpi(rowsFlag,'first')
                [preRowsTime, postRowsTime] = deal(5e-3,-4.8e-3);
                ROT = phantom_getRowOverTime(exper_dirCell_e, eventNumCell_e, CamMetaStruct, preRowsTime, postRowsTime);
            elseif strcmpi(rowsFlag,'last')
                [preRowsTime, postRowsTime] = deal(-4.8e-3, 5e-3);
                ROT = phantom_getRowOverTime(exper_dirCell_e, eventNumCell_e, CamMetaStruct, preRowsTime, postRowsTime);
            elseif strcmpi(rowsFlag,'all')
                [preRowsTime, postRowsTime] = deal(5e-3);
                ROT = phantom_getRowOverTime(exper_dirCell_e, eventNumCell_e, CamMetaStruct, preRowsTime, postRowsTime);
            end
        catch
            disp('');
        end
        first = ROT(1,:);
        last = ROT(end,:);
        newRows = [first(:)';last(:)'];
        rowStorage = cat(3,rowStorage,newRows);
    end
    
    
    %% correlate among first rows
    %--- correlate
    if strcmpi(rowsFlag,'first')
        allRelevantRows = permute(rowStorage(1,:,:),[3 2 1]);
    elseif strcmpi(rowsFlag,'last')
        allRelevantRows = permute(rowStorage(end,:,:),[3 2 1]);
    elseif strcmpi(rowsFlag,'all')
        %-- in this case it's the same as 'first'.
        allRelevantRows = permute(rowStorage(1,:,:),[3 2 1]);
    end
    lagsVec = zeros(N,0);
    dLagVec = zeros(N,0);
    RowsCropped = cell(N,1);
    spatialVecsCrop = cell(N,1);
    for i = 2:N
        preRow = allRelevantRows(i-1,:);
        postRow = allRelevantRows(i,:);
        
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
    allPointsX = [];
    EvIdx = [];
    EvntsLength = [];
    phediSlope = [];
    phediIntlLocPix = [];
    L = size(rowStorage,2);
    for i=1:N
        Xshifted = (1:L)+totalLagSUM(i);
        
        currentIdxs = PhediStructCellRelevant{i}.PhediData.phedisInitialLocsInPix;
        EvIdx = cat(1,EvIdx,ones(length(currentIdxs),1)*i);
        EvntsLength = cat(1,EvntsLength,length(currentIdxs));
        thisRoundSlope = PhediStructCellRelevant{i}.PhediData.slopeIncline';
        thisRoundSlope = round(mean(thisRoundSlope,2));     %-> a workaround for new data structures with the slope value for each frame.
        phediSlope = cat(1,phediSlope,thisRoundSlope);
        phediIntlLocPix = cat(1,phediIntlLocPix,PhediStructCellRelevant{i}.PhediData.phedisInitialLocsInPix);
        allPointsX = cat(1,allPointsX,Xshifted(currentIdxs)');
    end
    
    %% cluster by find peaks
    smthFac = 11;
    smthSignalFirst = smooth(allRelevantRows(1,:),smthFac);
    smthSignalLast = smooth(allRelevantRows(N,:),smthFac);
    % [pksMINfirst,locsMINfirst] = findpeaks(-smthSignalFirst);
    [pksMAXfirst,locsMAXfirst] = findpeaks(smthSignalFirst);
    % [pksMINlast,locsMINlast] = findpeaks(-smthSignalLast);
    [pksMAXlast,locsMAXlast] = findpeaks(smthSignalLast);
    
    % CloseMIN = interp1(locsMINlast,locsMINlast,locsMINfirst,'nearest');
    % meanValMIN = mean([CloseMIN(:),locsMINfirst(:)],2,'omitnan');
    
    CloseMAX = interp1(locsMAXlast,locsMAXlast,locsMAXfirst,'nearest');
    meanValMAX = mean([CloseMAX(:),locsMAXfirst(:)],2,'omitnan');
    
    %--- create bins for sorting phedis
    BinBorders = unique([0;meanValMAX(:);L+max(totalLagSUM)]);
    Steps = multipleStepFunction(BinBorders);
    XBins = Steps(allPointsX);
    
    
    phedigroupsIndx = XBins.*phediSlope;
    [~,~,phediGrouping] = unique(phedigroupsIndx);
    phediGroupingCell_beta = mat2cell(phediGrouping,EvntsLength',1);
    phediIntlLocByEv = mat2cell(phediIntlLocPix,EvntsLength',1);
    %% produce a slopes vector
    groupsANDslopes = [phediGrouping, phediSlope];
    [uniqeRows, ia, ~] = unique(groupsANDslopes,'rows');
    if size(uniqeRows,1)~=length(unique(phediGrouping))
        error('group number to slope sign sin''t a unique function');
    end
    slopesByOrder = uniqeRows(:,2);
    %     phediIntlLocPix(ia)
    
    %% check for duplicates in events and ignore
    for i=1:length(phediGroupingCell_beta)
        [occur,value] = hist(phediGroupingCell_beta{i},unique(phediGroupingCell_beta{i}));
        nonSingleVals = value(occur>1);
        for j=nonSingleVals(:)'
            phediGroupingCell_beta{i}(phediGroupingCell_beta{i}==j) = 0;
        end
    end
    phediGroupingCell = phediGroupingCell_beta;
    
    %% order of phedi groups
    phediGroupingMAT_beta = cell2mat(phediGroupingCell_beta);
    uniqueGroups = unique(phediGroupingMAT_beta);
    approxGroupLoc = zeros(length(uniqueGroups),1);
    for i=1:length(uniqueGroups)
        approxGroupLoc(i) = mean(allPointsX(phediGroupingMAT_beta==uniqueGroups(i)));
    end
    [~,I] = sort(approxGroupLoc);
    GroupsOrder = uniqueGroups(I);
    
end