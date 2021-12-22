%% plot phedi property by group
% Cell1 = load('2018-10-10\allPhediStructures_17186.mat');
% Cell2 = load('2018-10-10\allPhediStructures_165017.mat');
Cell1 = Cell1.PhediStructCell;
Cell2 = Cell2.PhediStructCell;
% load('E:\Frics\2018-10-10\allPhediStructures_17186.mat')
% PhediStructCellOfCell = {Cell1,Cell2,Cell3};
PhediStructCellOfCell = {Cell1,Cell2};
% PhediStructCellOfCell = {Cell2};
% PhediStructCellOfCell = {PhediStructCell};


%--- define damaged events
damaged_events20180829 = {...           % for 29-8-2018
    '12-8-52 15','12-8-52 17',...
    '14-2-11 6','14-2-11 12',...
    '14-45-14 10','14-45-14 19','14-45-14 28','14-45-14 12',...
    '15-9-25 19',...
    '15-28-23 29','15-28-23 19','15-28-23 18','15-28-23 23','15-28-23 13','15-28-23 8',...
    '18-1-9 14','18-1-9 11',...
    '18-20-0 24','18-20-0 19','18-20-0 7','18-20-0 4'...
    };
damaged_events_wDate20180829 = cellfun(@(A) ['2018-8-29 ',A],damaged_events20180829,'UniformOutput',0);

damaged_events20181010 = {...           % for 10-10-2018
    '12-2-22 28','12-2-22 20',...
    '14-28-36 4',...
    '14-48-1 12','14-48-1 21',...
    '20-28-52 14',...
    '20-48-22 16','20-48-22 3','20-48-22 5',...
    };
damaged_events_wDate20181010 = cellfun(@(A) ['2018-10-10 ',A],damaged_events20181010,'UniformOutput',0);

damaged_events = damaged_events_wDate20180829(:);
damaged_events(length(damaged_events_wDate20180829)+(1:length(damaged_events_wDate20181010))) = damaged_events_wDate20181010(:);

%--- keep only relevant events
PhediStructCellMega = [PhediStructCellOfCell{:}];
relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCellMega));
excludeNums = [];
for Ev = relevantEvents
    %-- eliminate rules:
    StrCheck = [PhediStructCellMega{Ev}.ExperimentData.ExpDate,' ',PhediStructCellMega{Ev}.ExperimentData.ExpHour,' ',num2str(PhediStructCellMega{Ev}.ExperimentData.eventNum)];
    damagedEv = nnz(strcmp(damaged_events,StrCheck))>0; % elimintae becaause damaged
    precurEv = phantom_isPrecursor(PhediStructCellMega{Ev}.BigPicRotStruct);   % elimintae becaause precoursor
    FastSlowCf = PhediStructCellMega{Ev}.PhediData.Cf<0 || PhediStructCellMega{Ev}.PhediData.Cf>1260;   % elimintae becaause Cf not in range
    if damagedEv || precurEv || FastSlowCf
        excludeNums = cat(1,excludeNums,Ev);
    end
end



FinalRelevantEv = setdiff(relevantEvents, excludeNums);
PhediStructCellMega_filter = PhediStructCellMega(FinalRelevantEv);



for iEv = 1:length(PhediStructCellMega_filter)
    PhediStructCellMega_filter{iEv}.PhediDataSG = PhediStructCellMega_filter{iEv}.PhediData;
    PhediStructCellMega_filter{iEv}.PhediData = PhediStructCellMega_filter{iEv}.PhediDataSimpSmth;
end


[phediGroupingCell, slopesByOrder, GroupsOrder, phediIntlLocByEv] = phedi_groupPhedisAcrossEvents_megaCell(PhediStructCellMega_filter);
GroupsOrder(GroupsOrder==0) = [];
% [rnkMaxCell, rnkP2pCell, maxVelCell, p2pVelCell] = phedi_VelocityComparison_megaCell(PhediStructCellOfCell);
[rnkMaxCell, rnkP2pCell, maxVelCell, p2pVelCell, stretchCell, rnkStrchCell, shiftsCell, rnkShiftCell,...
    GammaCell, LocalGammaCell, rnkLocalGCell, LocalKCell]...
    = phedi_VelocityComparison_megaCell(PhediStructCellMega_filter);


[Cd, Cs, Cr, nu, ro, E, mu, Gamma, PlaneStrain, tau_p, Xc0]=CrackSolutionMaterialProperties;
%% get common relevant events
%--- exclude precursurs from 'phediGroupingCell' as well
% occRnkMaxCell = ~cellfun(@isempty,rnkMaxCell);
% occPhediGroupingCell = ~cellfun(@isempty,phediGroupingCell);
% totalEventsLogical = occRnkMaxCell & occPhediGroupingCell;

phediGroupingCellFilter =  phediGroupingCell;
% phediGroupingCellFilter(~totalEventsLogical) = {[]};
%% create matrix of phedi groups (columns) and event (rows) where content will vary
phediGroupingMat = cell2mat(phediGroupingCellFilter(:));
Ngroups = max(phediGroupingMat(:));
MatBasis = nan(length(phediGroupingCellFilter),Ngroups);
%--- create vector for reordering
[As, Ai] = sort(GroupsOrder);
[Bs, Bi] = sort(sort(GroupsOrder));
reOrder = [];
reOrder(Ai) = Bi;
reOrderedSlopes = slopesByOrder(reOrder);

%% create velocity ranking mat:
if 0
    
    maxVelRnkMat = MatBasis;
    rnkMaxCellFilter =  rnkMaxCell;
    %     rnkMaxCellFilter(~totalEventsLogical) = {[]};
    for i=1:length(phediGroupingCellFilter)
        tempIndexes = phediGroupingCellFilter{i};
        tempInputs = rnkMaxCellFilter{i};
        tempInputs(tempIndexes==0) = [];
        tempIndexes(tempIndexes==0) = [];
        maxVelRnkMat(i,tempIndexes) = tempInputs;
    end
    maxVelRnkMatOrder = maxVelRnkMat(:,reOrder);
    
    figure;
    [h, hcb] = imagescWithNan(maxVelRnkMatOrder,jet,[0 0 0]);
    xlabel('phedi index'); ylabel('event number'); ylabel(hcb, 'vel rank');
    title({'velocity ranking',num2str(CellOfStruct{find(totalEventsLogical,1,'first')}.BigPicRotStruct.details(1:19))});
    
    
    %     figure;
    %     subplot(2,1,1);
    %     [h, hcb] = imagescWithNan(maxVelRnkMatOrder(:,reOrderedSlopes==-1),jet,[0 0 0]);
    %     title('-1'); xlabel('phedi index'); ylabel('event number'); ylabel(hcb, 'vel rank');
    %     subplot(2,1,2);
    %     [h, hcb] = imagescWithNan(maxVelRnkMatOrder(:,reOrderedSlopes==1),jet,[0 0 0]);
    %     title('+1'); xlabel('phedi index'); ylabel('event number'); ylabel(hcb, 'vel rank');
    %     suptitle({'velocity ranking',num2str(CellOfStruct{find(totalEventsLogical,1,'first')}.BigPicRotStruct.details(1:19))});
end

%% create velocity/Cf mat:
if 0
    maxVelCellFilter =  maxVelCell;
    %     maxVelCellFilter(~totalEventsLogical) = {[]};
    valVelMat = MatBasis;
    for i=1:length(phediGroupingCellFilter)
        tempIndexes = phediGroupingCellFilter{i};
        tempInputs = maxVelCellFilter{i};
        tempInputs(tempIndexes==0) = [];
        tempIndexes(tempIndexes==0) = [];
        valVelMat(i,tempIndexes) = tempInputs;
    end
    valVelMatOrder = valVelMat(:,reOrder);
    %--- normalize by event Cf
    %     CfByEvent = nan(size(valVelMatOrder,1),1);
    CfByEvent = cellfun(@(A) A.PhediData.Cf,PhediStructCellMega_filter);
    %--- normalize velocity
    valVelMatNormOrder = bsxfun(@rdivide,valVelMatOrder,CfByEvent(:));
    
    figure;
    [h, hcb] = imagescWithNan(valVelMatNormOrder,jet,[0 0 0]);
    title(['v_{max}/C_f ',num2str(CellOfStruct{find(totalEventsLogical,1,'first')}.BigPicRotStruct.details(1:19))]);
    xlabel('phedi index'); ylabel('event number'); ylabel(hcb, 'vel value');
    
    %     figure;
    %     subplot(2,1,1);
    %     [h, hcb] = imagescWithNan(valVelMatNormOrder(:,reOrderedSlopes==-1),jet,[0 0 0]);
    %     title('-1'); xlabel('phedi index'); ylabel('event number'); ylabel(hcb, 'vel value');
    %     subplot(2,1,2);
    %     [h, hcb] = imagescWithNan(valVelMatNormOrder(:,reOrderedSlopes==+1),jet,[0 0 0]);
    %     title('+1'); xlabel('phedi index'); ylabel('event number'); ylabel(hcb, 'vel value');
    %     suptitle({'velocity normalized by Cf',num2str(CellOfStruct{find(totalEventsLogical,1,'first')}.BigPicRotStruct.details(1:19))});
end


%% create p2p rank mat:
if 0
    rnkP2pCellFilter =  rnkP2pCell;
    %     rnkP2pCellFilter(~totalEventsLogical) = {[]};
    p2pRnkMat = MatBasis;
    for i=1:length(phediGroupingCellFilter)
        tempIndexes = phediGroupingCellFilter{i};
        tempInputs = rnkP2pCellFilter{i};
        tempInputs(tempIndexes==0) = [];
        tempIndexes(tempIndexes==0) = [];
        p2pRnkMat(i,tempIndexes) = tempInputs;
    end
    p2pRnkMatOrder = p2pRnkMat(:,reOrder);
    
    
    figure;
    [h, hcb] = imagescWithNan(p2pRnkMatOrder(:,reOrderedSlopes==-1),jet,[0 0 0]);
    xlabel('phedi index'); ylabel('event number'); ylabel(hcb, 'p2p ranking');
    title({'p2p ranking',num2str(CellOfStruct{find(totalEventsLogical,1,'first')}.BigPicRotStruct.details(1:19))});
    
    
    %     figure;
    %     subplot(2,1,1);
    %     [h, hcb] = imagescWithNan(p2pRnkMatOrder(:,reOrderedSlopes==-1),jet,[0 0 0]);
    %     title('-1'); xlabel('phedi index'); ylabel('event number'); ylabel(hcb, 'p2p ranking');
    %     subplot(2,1,2);
    %     [h, hcb] = imagescWithNan(p2pRnkMatOrder(:,reOrderedSlopes==1),jet,[0 0 0]);
    %     title('+1'); xlabel('phedi index'); ylabel('event number'); ylabel(hcb, 'p2p ranking');
    %     suptitle({'p2p ranking',num2str(CellOfStruct{find(totalEventsLogical,1,'first')}.BigPicRotStruct.details(1:19))});
end


%% by velocity
if 0
    figure;
    plot(rnkMaxMat,phediGroupingMat,'.');
    
    %-- velocity normalized by front velocity
    PhediStructCell_scrap = CellOfStruct(1:30);
    maxVelCell(cellfun(@isempty,maxVelCell))={0};
    cellLogical = cellfun(@(A) issame(A,0), maxVelCell);
    cellLogical = ~cellLogical;
    PhediStructCell_scrap = PhediStructCell_scrap(cellLogical);
    maxVelCell_scrap = maxVelCell(cellLogical);
    maxVelNormCell = cellfun(@(A,B) A./B.PhediData.Cf,maxVelCell_scrap,PhediStructCell_scrap,'UniformOutput',0);
    maxVelMat_scrap = cell2mat(maxVelCell_scrap);
    figure; plot(phediGroupingMat,maxVelMat_scrap,'.');
    maxVelNormMat= cell2mat(maxVelNormCell);
    figure; plot(phediGroupingMat,maxVelNormMat,'.');
end

%% plot by stretch factor
if 0
    stretchCellFilter = stretchCell;
    %     stretchCellFilter(~totalEventsLogical) = {[]};
    valStrchMat = MatBasis;
    for i=1:length(phediGroupingCellFilter)
        tempIndexes = phediGroupingCellFilter{i};
        tempInputs = stretchCellFilter{i};
        tempInputs(tempIndexes==0) = [];
        tempIndexes(tempIndexes==0) = [];
        valStrchMat(i,tempIndexes) = tempInputs;
    end
    valStrchMatOrder = valStrchMat(:,reOrder);
    figure; [h, hcb] = imagescWithNan(valStrchMatOrder,jet,[0 0 0]);
    xlabel('phedi index'); ylabel('event number'); ylabel(hcb, 'strtch factor');
    title('stretch factor');
    
end

%% create strtch factor rank mat:
if 0
    rnkStrchCellFilter =  rnkStrchCell;
    %     rnkStrchCellFilter(~totalEventsLogical) = {[]};
    rnkStrchMat = MatBasis;
    for i=1:length(phediGroupingCellFilter)
        tempIndexes = phediGroupingCellFilter{i};
        tempInputs = rnkStrchCellFilter{i};
        tempInputs(tempIndexes==0) = [];
        tempIndexes(tempIndexes==0) = [];
        rnkStrchMat(i,tempIndexes) = tempInputs;
    end
    rnkStrchMatOrder = rnkStrchMat(:,reOrder);
    figure;
    %     [h, hcb] = imagescWithNan(rnkStrchMatOrder,jet,[0 0 0]);
    rnkStrchMatOrder2 = rnkStrchMatOrder;
    rnkStrchMatOrder2(:,isnan(mean(rnkStrchMatOrder)))=[];
    h = imagesc(rnkStrchMatOrder2);
    colormap(jet);
    hcb = colorbar;
    xlabel('phedi index'); ylabel('event number'); ylabel(hcb, 'strtch fact. ranking');
    title('stretch factor ranking');
    
    [B,I] = sort(mean(rnkStrchMatOrder2));
    rnkStrchMatOrder_sorted = rnkStrchMatOrder2(:,I);
    figure; imagesc(rnkStrchMatOrder_sorted);
    
    colorbar;
    colormap(jet)
    title({'ranking of strtch factores','colmns sotrd by mean vals'})
    
    
    %     figure;
    %     subplot(2,1,1);
    %     [h, hcb] = imagescWithNan(rnkStrchMatOrder(:,reOrderedSlopes==-1),jet,[0 0 0]);
    %     title('-1'); xlabel('phedi index'); ylabel('event number'); ylabel(hcb, 'p2p ranking');
    %     subplot(2,1,2);
    %     [h, hcb] = imagescWithNan(rnkStrchMatOrder(:,reOrderedSlopes==1),jet,[0 0 0]);
    %     title('+1'); xlabel('phedi index'); ylabel('event number'); ylabel(hcb, 'p2p ranking');
    %     suptitle({'p2p ranking',num2str(CellOfStruct{find(totalEventsLogical,1,'first')}.BigPicRotStruct.details(1:19))});
end

%% rank local gamma
if 1
    rnkLocalGCellFilter =  rnkLocalGCell;
    %     rnkLocalGCellFilter(~totalEventsLogical) = {[]};
    rnkLclGMat = MatBasis;
    for i=1:length(phediGroupingCellFilter)
        tempIndexes = phediGroupingCellFilter{i};
        tempInputs = rnkLocalGCellFilter{i};
        tempInputs(tempIndexes==0) = [];
        tempIndexes(tempIndexes==0) = [];
        rnkLclGMat(i,tempIndexes) = tempInputs;
    end
    rnkLclGMatOrder = rnkLclGMat(:,reOrder);
    locG_mean = mean(rnkLclGMatOrder,'omitnan');
    locG_med = median(rnkLclGMatOrder,'omitnan');
    locG_std = std(rnkLclGMatOrder,'omitnan');
    [B,I] = sort(locG_mean);
    [B,I_med] = sort(locG_med);
    locG_mean_srtd = locG_mean(I);
    locG_med_srtd = locG_med(I);
    locG_std_srtd = locG_std(I);
    rnkLclGMatOrder_srtd = rnkLclGMatOrder(:,I);
    rnkLclGMatOrder_srtd_med = rnkLclGMatOrder(:,I_med);
    
    %--- index by location
    idxByLoc = 1:size(rnkLclGMatOrder,2);
    
    %--- plot rank by color
    Nphds = size(rnkLclGMatOrder,2);
    phediNum = 1:Nphds;
    eventNum = 1:size(rnkLclGMatOrder,1);
    figure;
    [h, hcb] = imagescWithNan(rnkLclGMatOrder_srtd,jet(Nphds),[0 0 0],phediNum([1,end]),eventNum([1,end]),1);
    g = gca;
    g.XTick = 1:length(I);
    g.XTickLabel = idxByLoc(I);
    xlabel('phedi index'); ylabel('event number'); ylabel(hcb, 'local G ranking');
    title('local \Gamma ranking sorted by mean');
    
    figure;
    [h, hcb] = imagescWithNan(rnkLclGMatOrder_srtd_med,jet(Nphds),[0 0 0],phediNum([1,end]),eventNum([1,end]),1);
    g = gca;
    g.XTick = 1:length(I_med);
    g.XTickLabel = idxByLoc(I_med);
    xlabel('phedi index'); ylabel('event number'); ylabel(hcb, 'local G ranking');
    title('local \Gamma ranking sorted by median');
    
    
    %--- plot mean vals and std
    figure;
    errorbar(phediNum,locG_mean_srtd, locG_std_srtd,'o','LineWidth',2);
    title('mean of local \Gamma ranking with std');
    
    %--- plot boxed of column distribution
    figure;
    boxplot(rnkLclGMatOrder_srtd_med, 'BoxStyle','filled');
    title('med of local \Gamma ranking');
    
    %     figure;
    %     subplot(2,1,1);
    %     [h, hcb] = imagescWithNan(rnkStrchMatOrder(:,reOrderedSlopes==-1),jet,[0 0 0]);
    %     title('-1'); xlabel('phedi index'); ylabel('event number'); ylabel(hcb, 'p2p ranking');
    %     subplot(2,1,2);
    %     [h, hcb] = imagescWithNan(rnkStrchMatOrder(:,reOrderedSlopes==1),jet,[0 0 0]);
    %     title('+1'); xlabel('phedi index'); ylabel('event number'); ylabel(hcb, 'p2p ranking');
    %     suptitle({'p2p ranking',num2str(CellOfStruct{find(totalEventsLogical,1,'first')}.BigPicRotStruct.details(1:19))});
end

%% plot by local G
if 0
    LocalGammaCellFilter = LocalGammaCell;
    %     LocalGammaCellFilter(~totalEventsLogical) = {[]};
    LocalGammaMat = MatBasis;
    for i=1:length(phediGroupingCellFilter)
        tempIndexes = phediGroupingCellFilter{i};
        tempInputs = LocalGammaCellFilter{i};
        tempInputs(tempIndexes==0) = [];
        tempIndexes(tempIndexes==0) = [];
        LocalGammaMat(i,tempIndexes) = tempInputs;
    end
    LocalGammaMatOrder = LocalGammaMat(:,reOrder);
    figure; [h, hcb] = imagescWithNan(real(LocalGammaMatOrder),jet,[0 0 0]);
    xlabel('phedi index'); ylabel('event number'); ylabel(hcb, 'local Gamma');
    title('local /Gamma');
    
    
    [B,I] = sort(mean(LocalGammaMatOrder,'omitnan'));
    LocalGammaMatOrder_sorted = LocalGammaMatOrder(:,I);
    figure; imagescWithNan(LocalGammaMatOrder_sorted,jet,[0 0 0]);
    title({'Local \Gamma','colmns sotrd by mean vals'})
    
    %--- plot mean vals and std
    locG_mean = mean(LocalGammaMatOrder_sorted,'omitnan');
    locG_std = std(LocalGammaMatOrder_sorted,'omitnan');
    figure;
    errorbar(1:(size(LocalGammaMatOrder_sorted,2)),locG_mean, locG_std,'o','LineWidth',2);
    title('mean of local \Gamma with srd');
    
    
end

%% plot by shift factor
if 0
    
    %     shiftsCell, rnkShiftCell
    
    shiftsCellFilter = shiftsCell;
    %     shiftsCellFilter(~totalEventsLogical) = {[]};
    valShftMat = MatBasis;
    for i=1:length(phediGroupingCellFilter)
        tempIndexes = phediGroupingCellFilter{i};
        tempInputs = stretchCellFilter{i};
        tempInputs(tempIndexes==0) = [];
        tempIndexes(tempIndexes==0) = [];
        valShftMat(i,tempIndexes) = tempInputs;
    end
    valShftMatOrder = valShftMat(:,reOrder);
    figure; [h, hcb] = imagescWithNan(valShftMatOrder,jet,[0 0 0]);
    xlabel('phedi index'); ylabel('event number'); ylabel(hcb, 'shift factor');
    title('shift factor');
    
end

%% create shift factor rank mat:
if 0
    rnkShftCellFilter =  rnkShiftCell;
    %     rnkStrchCellFilter(~totalEventsLogical) = {[]};
    rnkShftMat = MatBasis;
    for i=1:length(phediGroupingCellFilter)
        tempIndexes = phediGroupingCellFilter{i};
        tempInputs = rnkShftCellFilter{i};
        tempInputs(tempIndexes==0) = [];
        tempIndexes(tempIndexes==0) = [];
        rnkShftMat(i,tempIndexes) = tempInputs;
    end
    rnkShftMatOrder = rnkShftMat(:,reOrder);
    figure;
    [h, hcb] = imagescWithNan(rnkShftMatOrder,jet,[0 0 0]);
    xlabel('phedi index'); ylabel('event number'); ylabel(hcb, 'strtch fact. ranking');
    title('shift factor ranking');
    
    rnkShftMatOrder2 = rnkShftMatOrder;
    rnkShftMatOrder2(:,isnan(mean(rnkStrchMatOrder)))=[];
    [B,I] = sort(mean(rnkShftMatOrder2));
    rnkShftMatOrder_sorted = rnkShftMatOrder2(:,I);
    figure; imagesc(rnkShftMatOrder_sorted);
    
    colorbar;
    colormap(jet)
    title({'ranking of shift factores','colmns sotrd by mean vals'})
    
    %     figure;
    %     subplot(2,1,1);
    %     [h, hcb] = imagescWithNan(rnkStrchMatOrder(:,reOrderedSlopes==-1),jet,[0 0 0]);
    %     title('-1'); xlabel('phedi index'); ylabel('event number'); ylabel(hcb, 'p2p ranking');
    %     subplot(2,1,2);
    %     [h, hcb] = imagescWithNan(rnkStrchMatOrder(:,reOrderedSlopes==1),jet,[0 0 0]);
    %     title('+1'); xlabel('phedi index'); ylabel('event number'); ylabel(hcb, 'p2p ranking');
    %     suptitle({'p2p ranking',num2str(CellOfStruct{find(totalEventsLogical,1,'first')}.BigPicRotStruct.details(1:19))});
end

%% collapse phedi by stretch ranking
if 0
    %     phediGroupingCell, phediIntlLocByEv
    %     rnkStrchCell
    %     PhediStructCellMega_filter
    %     stretchCell
    
    %     [~,refEv] = min(cellfun(@(A) sum(A==0), phediGroupingCell));
    %     nnzInGrps = phediGroupingCell{refEv}(phediGroupingCell{refEv}~=0);
    %     refGroupNum = nnzInGrps(1);
    %     refPhediGroupNum = find(phediGroupingCell{refEv}==phediNum);
    
    %         ref_relativeGamma = LocalGammaCell{refEv}(:)';
    %         ref_groupNums = phediGroupingCell{refEv}(:)';
    
    
    
    
    phediGroupingMat = cell2mat(phediGroupingCellFilter);
    localGammamat = cell2mat(LocalGammaCell);
    ref_groupNums = unique(phediGroupingMat);
    ref_relativeGamma = nan(size(ref_groupNums));
        for u=1:length(ref_groupNums)
            ref_relativeGamma(u) = mean(localGammamat(ref_groupNums(u)==phediGroupingMat)); % assign a mean gamma for each phedi
    %         GammasTempVec = localGammamat(ref_groupNums(u)==phediGroupingMat);
    %         ref_relativeGamma(u) = GammasTempVec(1);
        end
    
%     for u=1:length(phediGroupingMat)
%         ref_relativeGamma4cell(u) = mean(localGammamat(phediGroupingMat(u)==phediGroupingMat)); % assign a mean gamma for each phedi
%     end
%     ref_relativeGamma4cell = ref_relativeGamma4cell(:);
%     ref_relativeCellByEv = mat2cell(ref_relativeGamma4cell,cellfun(@length,phediGroupingCellFilter),1);
    
    
    for i=1:length(PhediStructCellMega_filter)  % iterate on events
        PhediData = PhediStructCellMega_filter{i}.PhediData;
        PhediLocReduced = phedi_reduceSinFromLocation(PhediData);
        PhediLocReduced = bsxfun(@minus,PhediLocReduced,mean(PhediLocReduced(1:30,:),'omitnan'));
        spaceAxis = PhediData.x_mins_x_tip;
        inBlockLogicals = logical(PhediData.inBlockLogicals);
        Colors = MyVaryColor(size(PhediLocReduced,2));
        
        initialLocs_i = phediIntlLocByEv{i}(:)';
        groups_i = phediGroupingCell{i}(:)';
        
        
        
        %                 figure;
        %                 axes_colps = subplot(2,1,1);
        %                 hold on;
        %                 axes_reglr = subplot(2,1,2);
        %                 hold on;
        
        collectSapce_colps = [];
        collectPhediLoc_colps = [];
        collectSapce_reglr = [];
        collectPhediLoc_reglr = [];
        
        for p = 1:size(PhediLocReduced,2) % iterate on each phedi
            group_p = groups_i(PhediData.phedisInitialLocsInPix(p)==initialLocs_i);
            localGamma_p = ref_relativeGamma(ref_groupNums==group_p);   % take ref gamma
            if group_p==0 || isempty(localGamma_p)
                continue
            end
            
            v = PhediData.Cf;
            k=Cs/Cd;%Broberg p.330
            alpha_d=(1-(v./Cd).^2).^0.5;
            alpha_s=(1-(v./Cs).^2).^0.5;
            D=4*alpha_d.*alpha_s-(1+alpha_s.^2).^2;
            %----- Following Broberg p.334,336
            A=2*(1-k^2)*alpha_s.*v.^2./D/Cs^2;
            %             StrchBy = sqrt(A)./sqrt(mu*4*(1-k^2)./localGamma_p);    % generate stretching factor
            StrchBy = (sqrt(localGamma_p*mu*4*(1-k^2))/sqrt(A))*v*alpha_s/D; % generate stretching factor
            nonLocGStretch = (sqrt(GammaCell{i}(1)*mu*4*(1-k^2))/sqrt(A))*v*alpha_s/D;
            
            
            %             axes(axes_colps);
            %             plot(spaceAxis(inBlockLogicals(:,p),p),PhediLocReduced(inBlockLogicals(:,p),p)/sqrt(localGamma_p),...
            %                 '.-','Color',Colors(p,:),'LineWidth',1.5);
            %             axes(axes_reglr);
            %             plot(spaceAxis(inBlockLogicals(:,p),p),PhediLocReduced(inBlockLogicals(:,p),p),...
            %                 '.-','Color',Colors(p,:),'LineWidth',1.5);
            
            collectSapce_colps = cat(1,collectSapce_colps,spaceAxis(inBlockLogicals(:,p),p));
            collectPhediLoc_colps = cat(1,collectPhediLoc_colps,PhediLocReduced(inBlockLogicals(:,p),p)/sqrt(localGamma_p));
%                         collectPhediLoc_colps = cat(1,collectPhediLoc_colps,PhediLocReduced(inBlockLogicals(:,p),p)/StrchBy);
            collectSapce_reglr = cat(1,collectSapce_reglr,spaceAxis(inBlockLogicals(:,p),p));
            collectPhediLoc_reglr = cat(1,collectPhediLoc_reglr,PhediLocReduced(inBlockLogicals(:,p),p)/sqrt(GammaCell{i}(1)));
%                         collectPhediLoc_reglr = cat(1,collectPhediLoc_reglr,PhediLocReduced(inBlockLogicals(:,p),p)/nonLocGStretch);
            
            
        end
        
        dx = 1e-3;
        xshift = 5e-3;
        figure;
        %         axes_colps = subplot(2,1,1);
        %         yyaxis left
        plot(collectSapce_colps,collectPhediLoc_colps,...
            '.','Color','k','LineWidth',1.5,'DisplayName','colapsed data');
        hold on;
        plot(collectSapce_reglr,collectPhediLoc_reglr,...
            '.','Color',[1 1 1]*0.5,'LineWidth',1.5,'DisplayName','original data');
        yyaxis right
        [x_out_colps, var_out_colps, mean_out_colps] = my_movVar(collectSapce_colps,collectPhediLoc_colps,dx,xshift);
        %         plot(x_out_colps, mean_out_colps./sqrt(var_out_colps),'b','DisplayName','colapse \mu/\sigma');
        plot(x_out_colps, sqrt(var_out_colps)./mean_out_colps,'b','DisplayName','colapse \sigma/mean');
        [x_out_reglr, var_out_reglr, mean_out_reglr] = my_movVar(collectSapce_reglr,collectPhediLoc_reglr,dx,xshift);
        %         plot(x_out_reglr, mean_out_reglr./sqrt(var_out_reglr),'r','DisplayName','original \mu/\sigma');
        plot(x_out_reglr, sqrt(var_out_reglr)./mean_out_reglr,'r','DisplayName','original \sigma/mean');
        legend('show');
        title([PhediStructCellMega_filter{i}.BigPicRotStruct.details,' Cf=',num2str(PhediData.Cf)]);
        
        
        %         axes_reglr = subplot(2,1,2);
        %         yyaxis left
        %         plot(collectSapce_reglr,collectPhediLoc_reglr,...
        %             '.','Color','k','LineWidth',1.5);
        %         yyaxis right
        %         [x_out, var_out, mean_out] = my_movVar(collectSapce_reglr,collectPhediLoc_reglr,dx,xshift);
        %         plot(x_out, mean_out./sqrt(var_out),'b');
        
    end
    
end

%% mean of local gamma vs. gamma from sg
if 0
    locGammaMean = cellfun(@mean, LocalGammaCell);
    Gamma = cellfun(@mean, GammaCell);
    Cf_vec = cellfun(@(A) A.PhediData.Cf, PhediStructCellMega_filter);
    
    figure;
    scatter(Gamma,locGammaMean,80,Cf_vec,'o','filled','MarkerEdgeColor','r','DisplayName','mean(\Gamma_{local}) vs. \Gamma_{sg}');
    axis square
    hold on; plot([1 9],[1 9],'k','DisplayName','y=x');
    xlabel('\Gamma from sg measurement');
    ylabel('local \Gamma from phedi fitting');
    
    %-- find relation between Gamma and Cf. Cf_filter,G_filter are vecotrs
    %genetrated in the script 'ExtractDataCreateFigureFromDAandVstruct'
    GammaFit=slmengine(Cf_filter,G_filter,'plot','off','extrapolation','linear','concaveup','on','rightmaxslope',0,'rightminslope',0,'result','pp');
    GvsCf = ppval(GammaFit,Cf_vec);
    GammaTarget = 5;
    gammaRatio2_refCf = GvsCf/GammaTarget;
    scatter(Gamma,locGammaMean.*gammaRatio2_refCf,100,Cf_vec,'p','filled','MarkerEdgeColor',[0 0.8 0],'DisplayName','mean(\Gamma_{local}) fixed 2 \Gamma_{target by Cf} vs. \Gamma_{sg}');
    
    %--- get local effectiv Cf
    
    
    maxSrchRegn = 0.03;
    A_avgDist = 0.005;
    dxIntrp = 5e-5;
    intrpAx = -0.05:dxIntrp:0.05;
    maxVelOfPhdVec = cell(size(PhediStructCellMega_filter));
    CfVec = cell(size(PhediStructCellMega_filter));
    PeakVelOfMeanVec = cell(size(PhediStructCellMega_filter));
    
    for i=1:length(PhediStructCellMega_filter)
        PhediData = PhediStructCellMega_filter{i}.PhediData;
        spaceAxis4Vel = PhediData.x_mins_x_tip_4vel;
        phediVel= PhediData.PhediVelocity;
        Cf = PhediData.Cf;
        %--- get peak velocity of mean phedi
        Sall = phedi_FWHMofMean(PhediStructCellMega_filter{i},A_avgDist);
        VpeakXall = Sall.peakX;
        AvgPhediStructAll = Sall.AvgPhediStruct;
        [~,I] = min(abs(AvgPhediStructAll.X_mean- VpeakXall));
        peakVelFromMean = AvgPhediStructAll.Vel_mean_x(I);
        %--- get peak velocity of each phedi
        thisEvMaxVel = [];
        for p = 1:size(phediVel,2)
            nonInfNanLogic = ~(isinf(spaceAxis4Vel(:,p)) | isinf(phediVel(:,p)) | isnan(spaceAxis4Vel(:,p)) | isnan(phediVel(:,p)));
            intrpdVel = interp1(spaceAxis4Vel(nonInfNanLogic,p),phediVel(nonInfNanLogic,p),intrpAx);
            intrpdVelSmth = smooth(intrpdVel,round(A_avgDist/dxIntrp));
            maxVel = max(intrpdVelSmth(abs(intrpAx)<=maxSrchRegn));
            thisEvMaxVel = cat(1,thisEvMaxVel,maxVel);
        end
        
        %-- contain results in vectors
        maxVelOfPhdVec{i} = thisEvMaxVel;
        CfVec{i} = repmat(Cf,size(phediVel,2),1);
        PeakVelOfMeanVec{i} = repmat(peakVelFromMean,size(phediVel,2),1);
    end
    
    CfVsVpeak = empModel_CfVsVpeak_2();
    getEfectiveCf = @(X) ppval(CfVsVpeak,X);
    effectiveCf = cellfun(@(A) getEfectiveCf(A), maxVelOfPhdVec, 'UniformOutput',0);
    meanEfecCf = cellfun(@mean, effectiveCf);
    stdEfecCf = cellfun(@std, effectiveCf);
    Cf_real = cellfun(@mean, CfVec);
    effectiveGamma = ppval(GammaFit,meanEfecCf);
    
    %-- add effective gamma
    scatter(Gamma,effectiveGamma,80,Cf_vec,'s','filled','MarkerEdgeColor','b','DisplayName','\Gamma_{target by effectv Cf} vs. \Gamma_{sg}');
    gammaRatio2_refCf = GvsCf/effectiveGamma;
    scatter(Gamma,locGammaMean.*gammaRatio2_refCf,100,Cf_vec,'p','filled','MarkerEdgeColor','c','DisplayName','mean(\Gamma_{local}) fixed 2 \Gamma_{target by effectv Cf} vs. /Gamma_{sg}');
    
    %-- add efective gamma, each phedi
    EfecCfPerPhedi = cell2mat(effectiveCf(:));
    effctvGammaPerPhedi = ppval(GammaFit,EfecCfPerPhedi);
    GammaSgPerPhedi = cell2mat(GammaCell);
    GvsEfctGPerPhd = plot(GammaSgPerPhedi, effctvGammaPerPhedi,'o','Color',[1 1 1]*0.8,'DisplayName','local \Gamma_{effective by Cf} per phedi');
    uistack(GvsEfctGPerPhd,'bottom');
    
    %-- lengend
    LGND = legend('show');
    LGND.Location = 'northoutside';
    
    %-- mean effective Cf vs real Cf
    figure; errorbar(Cf_real,meanEfecCf,stdEfecCf,'o');
    xlabel('real Cf [m/s]'); ylabel('mean(effective Cf) [m/s]');
    title('mean effective Cf and std, calculated by mean of rach phedi using func. Vslip(Cf)')
    hold on; plot([0 1300],[0 1300],'k','DisplayName','y=x');
    
end

%% local dependance of gamma and Cf
if 0
    %-- purpose of this section is to create cell arrays of data devided into
    % phedi indexes instead of events.
    
    %     maxVel_cell2 = maxVelCell;
    maxVel_cell2 = maxVelOfPhdVec;
    
    LocalGammaCellFilter = LocalGammaCell;
    LocalKCellFilter = LocalKCell;
    maxVelCellFilter =  maxVel_cell2;
    valVelMat = MatBasis;
    LocalGammaMat = MatBasis;
    LocalKMat = MatBasis;
    gammaFromSg = MatBasis;
    CfPerPhedi_frmSg = MatBasis;
    for i=1:length(phediGroupingCellFilter)
        tempIndexes_vel = phediGroupingCellFilter{i};
        tempInputs_vel = maxVelCellFilter{i};
        tempInputs_vel(tempIndexes_vel==0) = [];
        tempIndexes_vel(tempIndexes_vel==0) = [];
        valVelMat(i,tempIndexes_vel) = tempInputs_vel;
        
        tempIndexes_G = phediGroupingCellFilter{i};
        tempInputs_G = LocalGammaCellFilter{i};
        tempInputs_G(tempIndexes_G==0) = [];
        tempIndexes_G(tempIndexes_G==0) = [];
        LocalGammaMat(i,tempIndexes_G) = tempInputs_G;
        
        tempIndexes_K = phediGroupingCellFilter{i};
        tempInputs_K = LocalKCellFilter{i};
        tempInputs_K(tempIndexes_K==0) = [];
        tempIndexes_K(tempIndexes_K==0) = [];
        LocalKMat(i,tempIndexes_K) = tempInputs_K;
        
        tempIndexes_Gsg = phediGroupingCellFilter{i};
        tempInputs_Gsg = GammaCell{i};
        tempInputs_Gsg(tempIndexes_Gsg==0) = [];
        tempIndexes_Gsg(tempIndexes_Gsg==0) = [];
        gammaFromSg(i,tempIndexes_Gsg) = tempInputs_Gsg;
        
        tempIndexes_Cf = phediGroupingCellFilter{i};
        tempInputs_Cf = CfVec{i};
        tempInputs_Cf(tempIndexes_Cf==0) = [];
        tempIndexes_Cf(tempIndexes_Cf==0) = [];
        CfPerPhedi_frmSg(i,tempIndexes_Cf) = tempInputs_Cf;
        
    end
    valVelMatOrder = valVelMat(:,reOrder);
    CfVsVpeak = empModel_CfVsVpeak_2();
    EfectiveCfMat= ppval(CfVsVpeak,valVelMatOrder);
    LocalGammaMatOrder = LocalGammaMat(:,reOrder);
    LocalKMatOrder = LocalKMat(:,reOrder);
    SgGammaMatOrder = gammaFromSg(:,reOrder);
    CfByEvent = cellfun(@(A) A.PhediData.Cf,PhediStructCellMega_filter);
    CfPerPhedi_frmSg_Order = CfPerPhedi_frmSg(:,reOrder);
    
    
    
    GammaFit = empModel_GammaVsCf();
    EfectiveGMat= ppval(GammaFit,EfectiveCfMat);
    
    figure; plot(SgGammaMatOrder,EfectiveGMat);
    xlabel('\Gamma from sg'); ylabel('\Gamma_{effective} from phedi pk vel');
    
    figure; plot(SgGammaMatOrder,LocalGammaMatOrder);
    xlabel('\Gamma from sg'); ylabel('\Gamma_{local} per phedi');
    
    figure; plot(CfByEvent,EfectiveGMat);   % these two parameters may be correlated
    xlabel('C_f from mesh'); ylabel('\Gamma_{effective} from phedi pk vel');
    
    figure; plot(CfByEvent,LocalGammaMatOrder);
    xlabel('C_f from mesh'); ylabel('\Gamma_{local}');
    
    GammaOfCf = ppval(GammaFit,CfPerPhedi_frmSg_Order);
    LocalGammaNrom = LocalGammaMatOrder./GammaOfCf;
    figure; plot(CfPerPhedi_frmSg_Order,LocalGammaNrom,'.-');   % these may also be correlated
    
end

%% last figures:
if 0
    K_vec = cellfun(@(A) A.solAtInter.K, PhediStructCellMega_filter);
    KLocalMean = cellfun(@mean, LocalKCell);
    KLocalVar = cellfun(@var, LocalKCell);
    figure; errorbar(K_vec,KLocalMean,KLocalVar,'o');
    hold on; plot([0 2.2e5],[0 2.2e5],'k');
    axis([0 2.5 0 2.5]*1e5);
    Cf_vec = cellfun(@(A) A.PhediData.Cf, PhediStructCellMega_filter);
    hold on;
    scatter(K_vec,KLocalMean,[],Cf_vec,'o','filled');
    CB = colorbar;
    xlabel('K_{II} from sg'); ylabel('K_{II} from phedis mean')
    
    figure; hold on;
    plot(Cf_vec,K_vec,'o','DisplayName','K from sg');
    plot(Cf_vec,KLocalMean,'o','DisplayName','K from phedis');
    
    
    [Cd, Cs, Cr, nu, ro, E, mu, Gamma, PlaneStrain, tau_p, Xc0]=CrackSolutionMaterialProperties;
    v = Cf_vec;
    k=Cs/Cd;%Broberg p.330
    alpha_d=(1-(v./Cd).^2).^0.5;
    alpha_s=(1-(v./Cs).^2).^0.5;
    D=4*alpha_d.*alpha_s-(1+alpha_s.^2).^2;
    %----- Following Broberg p.334,336
    A=2*(1-k^2)*alpha_s.*v.^2./D/Cs^2;
    G_sg = A.*(K_vec.^2)./(mu*4*(1-k^2));
    G_phedis = A.*(KLocalMean.^2)./(mu*4*(1-k^2));
    figure; hold on;
    plot(Cf_vec,G_sg,'o','DisplayName','G from sg');
    plot(Cf_vec,G_phedis,'o','DisplayName','G from phedis');
    figure; plot(G_sg,G_phedis,'o');
    axis square
    axis([0 12 0 12]);
    
end