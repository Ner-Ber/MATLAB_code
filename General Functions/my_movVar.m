% [x_out, var_out, mean_out] = my_movVar(x,v,dx,xshift)
%
% my_movVar will calculate the moving variance over a window of length dx.
% x doesn't nescesarily need to be evenly spaced
function [x_out, var_out, mean_out] = my_movVar(x,v,dx,xshift)
    LogicSelect = ~(isnan(x) | isinf(x));
    x_filt = x(LogicSelect);
    v_filt = v(LogicSelect);
    jumpingVec = (min(x_filt)-dx):xshift:(max(x_filt)+dx);
    x_out = nan(size(jumpingVec));
    var_out = nan(size(jumpingVec));
    mean_out = nan(size(jumpingVec));
    %-- omit infs:
    v_filt(abs(v_filt)==inf)=nan;
    for i=1:length(jumpingVec)
        thisWindow = [(jumpingVec(i)-dx),jumpingVec(i)];
        TheseLogic = x_filt>=thisWindow(1) & x_filt<=thisWindow(2);
        var_out(i) = var(v_filt(TheseLogic),'omitnan');
        mean_out(i) = mean(v_filt(TheseLogic),'omitnan');
        
        x_out(i) = mean(x_filt(TheseLogic));
    end
    
    
    
end