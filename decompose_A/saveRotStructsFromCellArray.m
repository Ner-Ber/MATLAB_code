function saveRotStructsFromCellArray(cellArray, folderPath)
    mkdir(folderPath);
    for ii=1:numel(cellArray)
        name = cellArray{ii}.BigPicRotStruct.details;
        name = regexprep(name, '\s*', '_');
        name = [strrep(name, '=', ''), '.txt'];
        saveRotDataToFile(cellArray{ii}.BigPicRotStruct, fullfile(folderPath, name));        
    end
end