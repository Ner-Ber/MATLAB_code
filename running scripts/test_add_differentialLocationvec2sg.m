%---*** update data structre so that time vectors of strain gages will be
%differential.
%Then recalculate the model stuff

for j=[6:11] % 2:11
    %% load the previous data structure
    loadName = ['allPhediStructures_Exp',num2str(j),'.mat'];
    load(loadName);
    %% update structures
    PhediStructCellUD = PhediStructCell;
    for i=1:30
        try
            DataStruct = PhediStructCell{i};
            if ~isfield(PhediStructCell{i},'theWorkspace')
                UpdatedDataStruct = sg_updateSpaceVectors(DataStruct);
                DataStructWSol = phedi_addLEFMandCohesive2Struct(UpdatedDataStruct);
                PhediStructCellUD{i} = DataStructWSol;
            end
        catch
        end
    end
    %% save result
    name=['allPhediStructures_Exp',num2str(j),'UD.mat'];
    save(name,'PhediStructCellUD','-v7.3');
    
end