function UpdatedDataStruct = sg_updateSpaceVectors(DataStruct,varargin)
    % UpdatedDataStruct = sg_updateSpaceVectors(DataStruct,CfSmooth)
    
    CfSmooth = setDefaults4function(varargin,10);
    
    UpdatedDataStruct = DataStruct;
    BigPicRotStruct = DataStruct.BigPicRotStruct;
    
    x_mins_x_tip_sg = zeros(size(DataStruct.SgData.t_mins_t_tips));
    for i = 1:size(DataStruct.SgData.t_mins_t_tips,2)
        %     x_min_Ttip_sgi = signal_ChangeTime2Space(...
        %         DataStruct.SgData.t_mins_t_tips(:,i), DataStruct.SgData.x_sg(i)*1e-3, BigPicRotStruct,CfSmooth);
        try
            x_min_Ttip_sgi = signal_ChangeTime2Space_mapping(...
                DataStruct.SgData.t_mins_t_tips(:,i), DataStruct.SgData.x_sg(i)*1e-3, BigPicRotStruct);
%             x_min_Ttip_sgi = signal_ChangeTime2Space_mapping20181912(...
%                             DataStruct.SgData.t_mins_t_tips(:,i), DataStruct.SgData.x_sg(i)*1e-3, BigPicRotStruct);
        catch
            disp();
        end
        x_mins_x_tip_sg(:,i) = x_min_Ttip_sgi(:);
    end
    %--- backup old and save new
    UpdatedDataStruct.SgData.x_mins_x_tips_constCf = UpdatedDataStruct.SgData.x_mins_x_tips;
    UpdatedDataStruct.SgData.x_mins_x_tips = x_mins_x_tip_sg;
    
end