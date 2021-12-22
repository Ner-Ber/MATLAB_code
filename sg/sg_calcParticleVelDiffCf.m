function particleVel = sg_calcParticleVelDiffCf(t_mins_t_tipOfSg, sgLoc, BigPicRotStruct, Uxx)
%x_min_Ttip = signal_ChangeTime2Space_mapping(t_mins_t_tipOfSg, sgLoc, BigPicRotStruct, Uxx)
%Uxx - strain in dL/L - NOT MILISTRAIN!!!
%sgLoc - IN METERS!
%
% signal_ChangeTime2Space_mapping will use signals t-t_tip and the front
% velocity to create vectors x-x_tip using a differential velocity.
% This will be done by using the time=f(location) function of the front as
% a lookup table. 


%% get mapping vectors
frontVel_interpMperS = BigPicRotStruct.frontVel_interpMperS;
t_front = movmean((BigPicRotStruct.frontTime_interp)/BigPicRotStruct.fps,2,'Endpoints','discard');
x_front = movmean(BigPicRotStruct.x,2,'Endpoints','discard');
%--- crop to keep only monotonic parts:
[~,locationIdx] = min(abs(x_front - sgLoc));
GradientTime = gradient(t_front);
bwlabelMonotonic = bwlabel(sign(GradientTime)==sign(GradientTime(locationIdx)));
monotonicRegion = bwlabelMonotonic == bwlabelMonotonic(locationIdx);
t_front_crop = t_front(monotonicRegion);
frontVel_crop = frontVel_interpMperS(monotonicRegion);
% x_front_crop = x_front(monotonicRegion);

%% find time of photoLocation:
[PhotoLocback,TtipsBigPic] = intersections(x_front,t_front,[sgLoc sgLoc],[-1e5 1e5]);

% particleVel = -interp1(t_front_crop-TtipsBigPic,x_front_crop-sgLoc,t_mins_t_tipOfSg,'linear');
Cf4partcleVel = interp1(t_front_crop-TtipsBigPic,frontVel_crop,t_mins_t_tipOfSg,'linear');
particleVel = -Cf4partcleVel.*Uxx;

