function [x_tip,t_tip] = phedi_find_t_tips_wPhotoLocation(BigPicRotStruct,PhotoLocation,spatialVec,...
                                    PhediLocation,measuredPhedisFromPlot,timeVec,varargin)
    % [x_tip,t_tip] = phedi_find_t_tips_wPhotoLocation(BigPicRotStruct,PhotoLocation,spatialVec,PhediLocation,measuredPhedisFromPlot,timeVec,Nsmth)
    %
    % phedi_find_t_tips_wPhotoLocation will find t_tips of the phedis given
    % the photograph location.
    
    %% set defaults and parameters:
    Nsmth = setDefaults4function(varargin,11);
    Nphd = length(measuredPhedisFromPlot);
    %%
    PhediLocOnBlock = bsxfun(@plus, PhediLocation,measuredPhedisFromPlot(:)')+(PhotoLocation-mean(spatialVec));
    %--- smooth phedi signals
    PhediLocOnBlockSmth = conv2(PhediLocOnBlock,ones(Nsmth,1)./Nsmth,'same');
    PhediLocOnBlockSmth([1:5,end-4:end],:) = PhediLocOnBlock([1:5,end-4:end],:);
    frontX = BigPicRotStruct.x;
    frontT = (BigPicRotStruct.frontTime_interp)/BigPicRotStruct.fps;
    t_tip = nan(Nphd,1);
    x_tip = nan(Nphd,1);
    for i = 1:Nphd
        [x_tip_i,t_tip_i] = intersections(PhediLocOnBlockSmth(:,i),timeVec,frontX,frontT);
        if length(t_tip_i)>1
            warning(['more than one intersection found in phedi num ',num2str(i)]);
            t_tip_i = mean(t_tip_i);
            x_tip_i = mean(x_tip_i);
        end
        t_tip(i) = t_tip_i;
        x_tip(i) = x_tip_i;
    end
    
    
end