function saveRotDataToFile(RotStruct, fileFullPath)
    fn = fieldnames(RotStruct);
    fileID = fopen(fileFullPath,'w');
    fprintf(fileID,'{\n');
    for kNum = 1:numel(fn)
        lastKey = kNum==numel(fn);
        k = fn{kNum};
        v = RotStruct.(k);
        fprintf(fileID,'\t\"%s\":', k);
        if ( isnumeric(v) )
            fprintf(fileID,'[\n');
            %--- split into main and last row for different notations
            v_main = v(1:end-1,:);
            v_end = v(end,:);
            
            %--- print main part, with commas after each line
            if ~isempty(v_main)
                sz_v = size(v_main);
                pString = repmat('%f, ', 1, sz_v(2));
                pString = pString(1:end-2);
                fprintf(fileID, ['\t\t[', pString, '],\n'], v_main');
            end
            %--- print last row, no comma at the end of the line
            sz_v = size(v_end);
            pString = repmat('%f, ', 1, sz_v(2));
            pString = pString(1:end-2);
            fprintf(fileID, ['\t\t[', pString, ']\n'], v_end');
            
            if lastKey
                fprintf(fileID,'\t]\n');
            else
                fprintf(fileID,'\t],\n');
            end
        else
            if lastKey
                fprintf(fileID,' \"%s\"\n', k);
            else
                fprintf(fileID,' \"%s\",\n', k);
            end
        end
        
    end
    fprintf(fileID,'\n}\n');
end