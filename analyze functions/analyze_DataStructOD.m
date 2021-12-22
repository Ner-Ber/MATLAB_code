function PhediStructCell = analyze_DataStructOD(ExperVec, Events)
% analyze_DataStructOD(ExperVec, EventVec)
%
% analyze_DataStructOD will build the DataStruct of the experiments and
% events inputed by usage of Movie_phedi_from_folder_2_data_defByName.
%
% ExperVec - vector of experiments in the current folder to analyze
% Events -  'all' - will analyze all events in the experiment
%           vector - will analyze the event specified in the vector for
%               each experiment
%           cell array - will contain either event vectors or the strinf
%               'all', each cell corresponds to a specific experiment by
%               its location.



% close all;
% clear;
exper = my_dir

%% iterate upon experiments
for ExperIdx=1:length(ExperVec)
    %--- prepare a cell array of events for reading
    if ~iscell(Events)
        if ischar(Events)
            AllFolders = dir_names([exper{ExperVec(ExperIdx)},'\acq132_094\multivent\']);
            eventsVec = 1:nnz(~cellfun(@isempty,strfind(AllFolders,'COOKED')));
        else
            eventsVec = Events;
        end
    else
        if ischar(Events{ExperIdx})
            AllFolders = dir_names([exper{ExperVec(ExperIdx)},'\acq132_094\multivent\']);
            eventsVec = 1:nnz(~cellfun(@isempty,strfind(AllFolders,'COOKED')));
        else
            eventsVec = Events{ExperIdx};
        end
    end
    
    %% itertate upon events
    PhediStructCell=cell(eventsVec(end),1);
    for EvNum = eventsVec(:)'
        try
            %--- don't run over existing measurements
            if length(PhediStructCell)>=EvNum
                if ~isempty(PhediStructCell{EvNum})
                    continue
                end
            end
            
            disp(['*** current event: ',num2str(EvNum),'***'])
            %--- get phedi data
            PhediStructCell{EvNum} = Movie_phedi_from_folder_2_data_defByName(ExperVec(ExperIdx),EvNum,...
                'preRowsTime',5e-3,'postRowsTime',5e-3,'prePhediTime',4e-3,'postPhediTime',4e-3,'QuickAndDirty',0);
            
        catch e
            fprintf(1,'The identifier was:\n%s',e.identifier);
            fprintf(1,'There was an error! The message was:\n%s',e.message);
        end
    end
    
    %% save
    name=['PhediStructures_',num2str(eventsVec(1)),'to',num2str(eventsVec(end)),'_Exp',num2str(ExperVec(ExperIdx)),'.mat'];
    save(name,'PhediStructCell','-v7.3');
end


end