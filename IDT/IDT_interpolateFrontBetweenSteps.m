function [frontIntepolatedStruct] = IDT_interpolateFrontBetweenSteps(frontStructure,varargin)
%[frontIntepolatedStruct] = IDT_interpolateFrontBetweenSteps(frontStructure)
% do better than finding step (by cubic interpolation)
% frontStructure - the outpus of the IDT_FindCfFromRowOverTime function.
% frontIntepolatedStruct - structure containing interpolated data.


[Nsmth] = setDefaults4function(varargin,7);

x = frontStructure.StepsPix;
t = frontStructure.StepTime;
xq = 1:length(frontStructure.fronLoc);
% tq = interp1(x,t,xq);
tq = interp1(x,t,xq,'pchip');
% tq = interp1(x,t,xq,'spline');

%--- smooth the front (assume slower changes)
tq_original = tq;
tq = smooth(tq,Nsmth);

dX=diff(xq);
dT=diff(tq);
FrontVel=dX(:)./dT(:);
VelLoc=movmean(xq,2,'Endpoints','discard');

%--- input to structure
frontIntepolatedStruct = struct;
frontIntepolatedStruct.xq = xq;
frontIntepolatedStruct.tq = tq;
frontIntepolatedStruct.tq_original = tq_original;
frontIntepolatedStruct.FrontVel = FrontVel;
frontIntepolatedStruct.VelLoc = VelLoc;