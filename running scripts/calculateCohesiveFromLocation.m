%%% plot location of phedis systematically
% close all;
% clear;
exper = my_dir

for ExperNum=1
    eventStart = 14;
    eventEnd = 23;
    
    moveOn = 1;
    failNum=0;
    EvNum=eventStart;
    
    PhediStructCell={};
    PhediLocsInterped={};
    PhediLocsInterpedZeroed={};
    minimalIdx={};
    cohesiveTime=[];
    cohesiveLength=[];
    PhediLocsInrp0dSumCell = {};
    
    while moveOn && EvNum<=eventEnd
        try
            if length(PhediStructCell)>=EvNum
                if ~isempty(PhediStructCell{EvNum})
                    EvNum = EvNum+1;
                    continue
                end
            end
            disp(['*** current event: ',num2str(EvNum),'***'])
            %--- get phedi data
            PhediStructCell{EvNum} = Movie_phedi_from_folder_2_data_defByName(ExperNum,EvNum,...
                'preRowsTime',3e-3,'postRowsTime',3e-3,'prePhediTime',4e-3,'postPhediTime',4e-3,'QuickAndDirty',0);
            
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
    
    %% save
    name=['PhediStructures_14_2_23_',num2str(ExperNum),'.mat'];
    save(name,'PhediStruct');
end

%% calculate MAYBE cohesive zone
% isEmpty = cellfun(@isempty,PhediStruct);
% EventNums = find(~isEmpty);
% PhediStructMinimal = PhediStruct(~isEmpty);
% t_mins_t_tip_Minimal = cellfun(@(x) bsxfun(@minus,x.PhediData.timeVec,x.PhediData.t_tips(:)'),PhediStructMinimal,'UniformOutput',0);
% t_mins_t_tip_sortdMinimal = cellfun(@(x) sort(x(:)),t_mins_t_tip_Minimal,'UniformOutput',0);
%
% PhediLocsInterped=cell(size(PhediStructMinimal));
% PhediLocsInterpedZeroed=cell(size(PhediStructMinimal));
% for cellIdx=1:length(PhediStructMinimal)
%     N_phedis = size(PhediStructMinimal{cellIdx}.PhediData.PhediLocation,2);
%     PhediLocsInterped{cellIdx} = nan(length(t_mins_t_tip_sortdMinimal{cellIdx}),N_phedis);
%     for i=1:N_phedis
%         PhediLocsInterped{cellIdx}(:,i)=interp1(t_mins_t_tip_Minimal{cellIdx}(:,i),...
%             PhediStructMinimal{cellIdx}.PhediData.PhediLocation(:,i),...
%             t_mins_t_tip_sortdMinimal{cellIdx},...
%             'linear','extrap');
%     end
%     PhediLocsInterpedZeroed{cellIdx} = bsxfun(@minus,PhediLocsInterped{cellIdx},mean(PhediLocsInterped{cellIdx}(300:350,:),1,'omitnan'));
% end
% PhediLocsInterpedStdMinimal = cellfun(@(x) std(x,0,2,'omitnan') , PhediLocsInterpedZeroed,'UniformOutput',0);
% PhediLocsInterpedStdMinimalSmooth250 = cellfun(@(x) smooth(x,250), PhediLocsInterpedStdMinimal,'UniformOutput',0);
% PhediLocsInterpedStdNormalized = cellfun(@(y,x) y./find_y0_at_x0(0,y,x), PhediLocsInterpedStdMinimalSmooth250,t_mins_t_tip_sortdMinimal,'UniformOutput',0);
% typicalTimes = cellfun(@(y,x) find_x0_at_y0(1/exp(1),y(x<=1e-4),x(x<=1e-4)), PhediLocsInterpedStdNormalized,t_mins_t_tip_sortdMinimal);
% Cf_vec = cellfun(@(x) x.PhediData.Cf,PhediStructMinimal);
% typicalLength = typicalTimes.*Cf_vec;
% figure; plot(Cf_vec/1255,abs(typicalLength),'r*');
%
%
