function PhediStruct = analyze_batchAnalyzePhedis(experVec,EventStart)

for ExperNum=experVec
    moveOn = 1;
    failNum=0;
    if nargin<2
        EvNum=1;
    else
        EvNum=EventStart;
    end
    
    PhediStruct={};
    
    while moveOn && EvNum<=45
        try
            if length(PhediStruct)>=EvNum
                if ~isempty(PhediStruct{EvNum})
                    EvNum = EvNum+1;
                    continue
                end
            end
            disp(['*** current event: ',num2str(EvNum),'***'])
            %--- get phedi data
%             PhediStruct{EvNum} = Movie_phedi_from_folder_2_data(ExperNum,EvNum,[],[],[],5e-3,5e-3,5e-3,5e-3);
%             PhediStruct{EvNum} = Movie_phedi_from_folder_2_data_defByName(ExperNum,EvNum,...
%                 'preRowsTime',5e-3,'postRowsTime',5e-3,'prePhediTime',3e-3,'postPhediTime',3e-3);
            PhediStruct{EvNum} = phedi_createStructureWAdd(ExperNum,EvNum,...
                'preRowsTime',5e-3,'postRowsTime',5e-3,'prePhediTime',3e-3,'postPhediTime',1.8e-3);
            
            EvNum = EvNum+1;
            failNum = 0;
        catch e
            fprintf(1,'The identifier was:\n%s',e.identifier);
            fprintf(1,'There was an error! The message was:\n%s',e.message);
            EvNum = EvNum+1;
            failNum=failNum+1;
            if failNum>=5
                break
            end
        end
        
    end
end