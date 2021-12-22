%% purpose
% the purposeof this script is to load and compare data from different
% lines in the same event.

%% load
matObj = matfile('allPhediStructures_singleLines_17-18-6.mat');
% matObj = matfile('allPhediStructures_specLines_sepLines_17-18-6');
% multipLineStruct = matObj.PhediStructCell(1,1);

%% create 3D ROT for events
names = fieldnames(matObj.PhediStructCell(1,1));
DataMat3DCell = cell(1,30);
for i=1:30; DataMat3DCell{i}=[]; end
for i = 1:numel(names)
    for c = find(~cellfun(@isempty,matObj.PhediStructCell(1,1).(names{i})))
        try
            DataMat3DCell{c} = cat(3,DataMat3DCell{c},matObj.PhediStructCell(1,1).(names{i}){c}.BigPicRotStruct.DataMatNorm);
        catch
            try
                DataMat3DCell{c} = cat(3,DataMat3DCell{c},matObj.PhediStructCell(1,1).(names{i}){c}.theWorkspace.BigPicRotStruct.DataMatNorm);
            catch
            end
        end
    end
end