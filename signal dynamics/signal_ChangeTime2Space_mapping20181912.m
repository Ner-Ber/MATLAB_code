function [x_min_Ttip, relevantLocations] = signal_ChangeTime2Space_mapping20181912(t_mins_t_tipOfSignal, PhotoLocation, BigPicRotStruct)
    %x_min_Ttip = signal_ChangeTime2Space_mapping(t_mins_t_tipOfSignal, PhotoLocation, BigPicRotStruct)
    %
    % signal_ChangeTime2Space_mapping will use signals t-t_tip and the front
    % velocity to create vectors x-x_tip using a differential velocity.
    % This will be done by using the time=f(location) function of the front as
    % a lookup table.
    
    
    %% get mapping vectors
    t_front = (BigPicRotStruct.frontTime_interp)/BigPicRotStruct.fps;
    x_front = BigPicRotStruct.x;
    %--- crop to keep only monotonic parts:
    [~,locationIdx] = min(abs(x_front - PhotoLocation));
    GradientTime = gradient(t_front);
    bwlabelMonotonic = bwlabel(sign(GradientTime)==sign(GradientTime(locationIdx)));
    monotonicRegion = bwlabelMonotonic == bwlabelMonotonic(locationIdx);
    t_front_crop = t_front(monotonicRegion);
    x_front_crop = x_front(monotonicRegion);
    
    %% find time of photoLocation:
    [PhotoLocback,TtipsBigPic] = intersections(x_front,t_front,[PhotoLocation PhotoLocation],[-1e5 1e5]);
    
    TT = t_front_crop-TtipsBigPic;
    LL = logical(TT);
    if ~~nnz(LL)        % the case where there are different time values
        x_min_Ttip = -interp1(t_front_crop-TtipsBigPic,x_front_crop-PhotoLocation,t_mins_t_tipOfSignal,'linear','extrap');
%         x_min_Ttip = interp1(t_front_crop-TtipsBigPic,(x_front_crop-PhotoLocation),-t_mins_t_tipOfSignal,'linear','extrap');
%         x_min_Ttip = interp1(-(t_front_crop-TtipsBigPic),(x_front_crop-PhotoLocation),-t_mins_t_tipOfSignal,'linear','extrap');
        
    else                % if the whole time vector is the same,   THIS IS A CHEATING BYPASS, NOT ACCURATE AT ALL
        XX = x_front_crop-PhotoLocation;
        if min(t_mins_t_tipOfSignal)>max(t_front_crop-TtipsBigPic)
            x_min_Ttip = XX(end);
        else
            x_min_Ttip = XX(1);
        end
    end
    
    relevantLocations = t_mins_t_tipOfSignal<max(t_front_crop-TtipsBigPic) & t_mins_t_tipOfSignal>min(t_front_crop-TtipsBigPic);
end