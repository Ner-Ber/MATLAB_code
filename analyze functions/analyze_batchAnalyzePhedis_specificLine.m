function PhediStructMultipLines = analyze_batchAnalyzePhedis_specificLine(ExperNum, events, lNum_cell)


PhediStructMultipLines = struct;
for l = 1:length(lNum_cell)
    lNum = lNum_cell{l};
    failNum=0;

    PhediStruct={};

    for EvNum=events
        try
%             if length(PhediStruct)>=EvNum
%                 if ~isempty(PhediStruct{EvNum})
%                     EvNum = EvNum+1;
%                     continue
%                 end
%             end
            disp(['*** current event: ',num2str(EvNum),'***'])
            %--- get phedi data
%             PhediStruct{EvNum} = Movie_phedi_from_folder_2_data(ExperNum,EvNum,[],[],[],5e-3,5e-3,5e-3,5e-3);
%             PhediStruct{EvNum} = Movie_phedi_from_folder_2_data_defByName(ExperNum,EvNum,...
%                 'preRowsTime',5e-3,'postRowsTime',5e-3,'prePhediTime',3e-3,'postPhediTime',3e-3);
            PhediStruct{EvNum} = phedi_createStructureWAdd(ExperNum,EvNum,...
                'preRowsTime',5e-3,'postRowsTime',5e-3,'prePhediTime',3e-3,'postPhediTime',1.8e-3, 'lineNum', lNum);

            failNum = 0;
        catch e
            fprintf(1,'The identifier was:\n%s\n',e.identifier);
            fprintf(1,'There was an error! The message was:\n%s\n',e.message);
            failNum=failNum+1;
            fprintf(1,'fail %d in a row\n',failNum);
        end

    end
    PhediStructMultipLines.(['line_',strrep(num2str(lNum),' ','_')]) = PhediStruct;
end