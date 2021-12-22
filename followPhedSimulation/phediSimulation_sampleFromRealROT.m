function [RowOverTime_preDamp, sampleLocations] = phediSimulation_sampleFromRealROT(...
        RowOverTime_true, x_true, frameLength, simpleSamp, FillFac)
    
    
    sampleLocations = ismember(x_true,1:frameLength); %~mod(x_true,1); % find pixle locations
    if simpleSamp   % this case is just undersampling the data on the x direction
        RowOverTime_preDamp = RowOverTime_true(:,sampleLocations);
    else    % this option will integrate by a given fill factor
        x_true_mod = mod(x_true,1);
        D = abs(sqrt(FillFac))/2;
        fillLogic = x_true_mod<=D | x_true_mod>=(1-D);
        RowOverTime_fill = RowOverTime_true(:,fillLogic);
        ROT_reshape = reshape(RowOverTime_fill,size(RowOverTime_fill,1),[],frameLength);
        RowOverTime_preDamp = permute(sum(ROT_reshape,2),[1,3,2])./(size(ROT_reshape,2));
    end
    
    
    
end