% relevantCells = ~cellfun(@(A) isfield(A,'theWorkspace'),PhediStructCell);
% PhediStructCellUpdated = PhediStructCell;
PhediStructCellUpdated = cell(size(PhediStructCell));
relevantCells = zeros(size(PhediStructCell));
% for i = find(relevantCells)
for i = 1:length(PhediStructCellUpdated)
    if isfield(PhediStructCell{i},'theWorkspace');
        try
            PhediStructCellUpdated{i} = sg_updateSpaceVectors(PhediStructCell{i}.theWorkspace.DataStructWPhedis);
        catch
            PhediStructCellUpdated{i} = PhediStructCell{i}.theWorkspace.DataStructWPhedis;
        end
    elseif phantom_isPrecursor(PhediStructCell{i}.BigPicRotStruct)
        PhediStructCellUpdated{i} = PhediStructCell{i};
    else
        exp_dir = PhediStructCell{i}.ExperimentData.ExpHour;
        Param = PhediStructCell{i}.Param;
        Param.frontFallDetermn = 2;
        BigPicRowOverTimeStruct = phantom_BigPicStruct(exp_dir, i,Param.preRowsTime,1e-10,Param.postRowsTime,1,1,'all',Param.frontFallDetermn);
        DataStruct.Param = Param;
        DataStruct.CamMeta = PhediStructCell{i}.CamMeta;
        DataStruct.ExperimentData = PhediStructCell{i}.ExperimentData;
        DataStruct.AsperityData = PhediStructCell{i}.AsperityData;
        DataStruct.BigPicRotStruct = BigPicRowOverTimeStruct;
        DataStruct.PhediData = PhediStructCell{i}.PhediDataPreSmth;
        rupture_kind = 'slow';
        DataStruct = sg_getDataAndFixShearSens(DataStruct,rupture_kind);
        %--- update phedi location and phptoLocation
        [bigPicLocationPix, ~] = Movie_syncLocationOf2Cameras(BigPicRowOverTimeStruct.t,BigPicRowOverTimeStruct.DataMatNorm,...
                DataStruct.AsperityData.timeCount,DataStruct.AsperityData.RowOverTime,[],[],1);
        PhotoLocation = BigPicRowOverTimeStruct.x(bigPicLocationPix);
        [x_tips,t_tips] = phedi_find_t_tips_wPhotoLocation(BigPicRowOverTimeStruct,PhotoLocation,DataStruct.AsperityData.spatialVec,...
                                DataStruct.PhediData.PhediLocation,DataStruct.PhediData.measuredPhedisFromPlot,DataStruct.PhediData.timeVec);
        t_mins_t_tip = bsxfun(@minus,DataStruct.PhediData.timeVec(:),t_tips(:)');
        DataStruct.PhediData.PhotoLocation = PhotoLocation;
        DataStruct.PhediData.t_tips = t_tips;
        DataStruct.PhediData.x_tips = x_tips;
        DataStruct.PhediData.t_mins_t_tip = t_mins_t_tip;
        
        DataStructWLocation = phedi_add_phediLocationVectors(DataStruct);
        %--- add velocity vectors:
        DataStructWVelocity = phedi_add_velocity_to_struct(DataStructWLocation);
        %--- smooth functions measured:
        DataStructWSmooth = phedi_addSmoothedFunctions2Struct(DataStructWVelocity);
        %--- add:
        DataStructWSol = phedi_addLEFMandCohesive2Struct(DataStructWSmooth,Param);
        
        
        PhediStructCellUpdated{i} = DataStructWSol;
        
%         tmpStruct = rmfield(PhediStructCell{i},{'CohesiveModelStruct','solAtInter','solAtSG'});
%         tmpPhediStruct = 
%         PhediStructCellUpdated{i} = phedi_addLEFMandCohesive2Struct(tmpStruct,tmpStruct.Param);
        relevantCells(i) = 1;
    end
end
relevantCells = ~~relevantCells;

% for i = find(relevantCells)
%     figure; hold on;
%     plot(PhediStructCellUpdated{i}.SgData.x_mins_x_tips(:,8),1e-3*(PhediStructCellUpdated{i}.SgData.Uxx(:,8)-mean(PhediStructCellUpdated{i}.SgData.Uxx(1:100,8))),'.-','Color',rgb('Pink'));
%     plot(PhediStructCellUpdated{i}.solAtSG.x,PhediStructCellUpdated{i}.solAtSG.Uxx,'r');
%     plot(PhediStructCell{i}.solAtSG.x,PhediStructCell{i}.solAtSG.Uxx,'b');
%     title({PhediStructCell{i}.BigPicRotStruct.details,...
%         ['\Gamma_{new}=',num2str(PhediStructCellUpdated{i}.solAtSG.Gamma),'  \Gamma_{old}=',num2str(PhediStructCell{i}.solAtSG.Gamma)]});
%     xlim([-0.05 0.05])
% end
%
%
