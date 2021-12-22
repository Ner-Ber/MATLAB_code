function [x_min_Ttip,relevantLocations] = signal_ChangeTime2Space(t_mins_t_tipOfSignal, PhotoLocation, BigPicRotStruct,varargin)
    % x_min_Ttip = signal_ChangeTime2Space(t_mins_t_tipOfSignal, PhotoLocation, BigPicRotStruct,smoothWidth,smoothType)
    %
    % signal_ChangeTime2Space will use signals t-t_tip and the front velocity
    % to create vectors x-x_tip. This will be done by using
    % f[x,t]=f[x-int(Cf*dt)].
    %
    % OPTIONAL:
    % smoothWidth = 'gaussian' (default), 'rect', 'none'
    % smoothType =  a scalar. 0-none; 15-default
    %
    
    
    
    %% set defaults:
    [smoothWidth,smoothType] = setDefaults4function(varargin,10,'gaussian');
    
    %% load data
    xFront = BigPicRotStruct.x;
    tFront = (BigPicRotStruct.frontTime_interp)/BigPicRotStruct.fps;
    frontVel_interpMperS = BigPicRotStruct.frontVel_interpMperS;
    frontVelLoc_interpM = BigPicRotStruct.frontVelLoc_interpM;
    
    
    %% prepare relevant data
    %--- time vector for velocity of front
    frontVelTime_interpS = interp1(xFront(:),tFront(:),frontVelLoc_interpM);
    %-- t_tip from Big Picture
    [PhotoLocback,TtipsBigPic] = intersections(xFront,tFront,[PhotoLocation PhotoLocation],[-1e5 1e5]);
    
    %% monotinoc region around photo location
    [~,locationIdx] = min(abs(frontVelLoc_interpM - PhotoLocation));
    GradientVel = gradient(frontVelTime_interpS);
    bwlabelMonotonic = bwlabel(sign(GradientVel)==sign(GradientVel(locationIdx)));
    monotonicRegion = bwlabelMonotonic == bwlabelMonotonic(locationIdx);
    
    %% time vector for each phedi
    frontCfTime_min_t_tip = frontVelTime_interpS(monotonicRegion) - TtipsBigPic;
    relevantCf = frontVel_interpMperS(monotonicRegion);
    Cf_for_timeVec = interp1(frontCfTime_min_t_tip,relevantCf,t_mins_t_tipOfSignal,'linear');
    %--- replace NaNs
    relevantLocations = ~isnan(Cf_for_timeVec);
    NanGroups = bwlabel(isnan(Cf_for_timeVec));
    Cf_for_timeVec(NanGroups==1) = relevantCf(1);
    Cf_for_timeVec(NanGroups==max(NanGroups)) = relevantCf(end);
    
    
    %% smooth Cf function
    
    if smoothWidth==0 || strcmpi(smoothType,'none')
        CfSmoothed  = Cf_for_timeVec;
    elseif strcmpi(smoothType,'gaussian')
        N = 4*round(abs(smoothWidth));
        kernel = pdf('Normal',-N:N,0,smoothWidth);
        CfSmoothed = conv(Cf_for_timeVec,kernel,'same');
    elseif strcmpi(smoothType,'rect')
        CfSmoothed = smooth(Cf_for_timeVec,smoothWidth);
    end
    
    %% create spatial vector by integration
    
    dt = mean(gradient(t_mins_t_tipOfSignal));
    x_min_Ttip_Shifted = -cumsum(dt.*CfSmoothed);
    
    %% constrain zero to stay put
    %--- when computing 'x_min_Ttip_Shifted' you get an inherentic shift becase
    %of differentiating then integrating.
    if ~(nnz(isnan(t_mins_t_tipOfSignal))==length(t_mins_t_tipOfSignal) || nnz(isnan(x_min_Ttip_Shifted))==length(x_min_Ttip_Shifted))
        dx = interp1(t_mins_t_tipOfSignal,x_min_Ttip_Shifted,0);
        x_min_Ttip = x_min_Ttip_Shifted-dx;
    else
        x_min_Ttip = x_min_Ttip_Shifted;
    end
    
    
end