%% plot phedi property by group
% Cell1 = load('2018-10-10\allPhediStructures_17186.mat');
% Cell2 = load('2018-10-10\allPhediStructures_165017.mat');
% Cell1 = Cell1.PhediStructCell;
% Cell2 = Cell2.PhediStructCell;
% load('E:\Frics\2018-10-10\allPhediStructures_17186.mat')
% PhediStructCellOfCell = {Cell1,Cell2,Cell3};
% PhediStructCellOfCell = {Cell1,Cell2};
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
% PhediStructCellMega = [PhediStructCellOfCell{:}];
PhediStructCellMega = PhediStructCell;
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

