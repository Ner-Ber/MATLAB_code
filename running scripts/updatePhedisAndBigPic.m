%% update BigPicStruct and followinf phedis


D = dir_names;
Experiments = find(~cellfun(@isempty,strfind(D,'_')));
for i=Experiments
    %% load data
    load(D{i});
    relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
    relevantEvents = relevantEvents(relevantEvents~=1);
    ThisExp = PhediStructCell{relevantEvents(1)}.ExperimentData.ExpHour;
    PhediStructCellOld = PhediStructCell;
    %% iterate events
    for iEv = relevantEvents
        %% create new BigPic
        DataStruct = PhediStructCellOld{iEv};
        Param = PhediStructCellOld{iEv}.Param;
        BigPicRotStruct = phantom_BigPicStruct(ThisExp,iEv,5e-3,[],5e-3,[],[],[],5);
        DataStruct.BigPicRotStruct = BigPicRotStruct;
        %% update stuff
        NewDataStruct = phedi_addOrChngPhediTtips(DataStruct);
        rupture_kind = 'slow';
        DataStructWsg = sg_getDataAndFixShearSens(NewDataStruct,rupture_kind);
        %--- add location vectors:
        DataStructWLocation = phedi_add_phediLocationVectors(DataStructWsg);
        %--- add velocity vectors:
        DataStructWVelocity = phedi_add_velocity_to_struct(DataStructWLocation);
        %--- smooth functions measured:
        DataStructWSmooth = phedi_addSmoothedFunctions2Struct(DataStructWVelocity);
        %-- ad smoothed by Savitzky Golay:
        DataStructSgSmooth = phedi_smoothWsgolay(DataStructWSmooth);
        %--- add:
        DataStructWSol = phedi_addLEFMandCohesive2Struct(DataStructSgSmooth,Param);
        
        PhediStructCell{iEv} = DataStructWSol;
    end
    %% resave
    hourStr = exp_getHourFromDataStructCell(PhediStructCell,i);
    name=['allPhediStructures_',hourStr,'.mat'];
    save(name,'PhediStructCell','-v7.3');
end