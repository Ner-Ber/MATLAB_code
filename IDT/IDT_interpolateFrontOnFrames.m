function [frontIntepolatedStruct] = IDT_interpolateFrontOnFrames(frontStructure)
%[frontIntepolatedStruct] = IDT_interpolateFrontBetweenSteps(frontStructure)
% do better than finding step (by cubic interpolation)
% frontStructure - the outpus of the IDT_FindCfFromRowOverTime function.
% frontIntepolatedStruct - structure containing interpolated data.

x = frontStructure.StepsPix;
t = frontStructure.StepTime;
%% find monotonic part of front
st_diff = sign(diff(t));    % give sign for increasing or decreasing region
i = find(diff(st_diff));    % find where regions switch
n = [i(:)' numel(st_diff)] - [0 i(:)']; % how long is each region
[~, LongIdx] = max(n);            % longest region 
I = [0 i(:)']+1;                % starting idexes of regions
monotonRegion = I(LongIdx):(I(LongIdx+1)-1);
%% crop for interpolation
x_crop = x(monotonRegion);
t_crop = t(monotonRegion);


%% interpolate
tq = t_crop(1):t_crop(end);
% xq = interp1(t_crop,x_crop,tq);
xq = interp1(t_crop,x_crop,tq,'pchip');
% xq = interp1(t_crop,x_crop,tq,'spline');

dX=diff(xq);
dT=diff(tq);
FrontVel=dX./dT;
VelLoc=movmean(xq,2,'Endpoints','discard');

%--- input to structure
frontIntepolatedStruct = struct;
frontIntepolatedStruct.xq = xq;
frontIntepolatedStruct.tq = tq;
frontIntepolatedStruct.FrontVel = FrontVel;
frontIntepolatedStruct.VelLoc = VelLoc;