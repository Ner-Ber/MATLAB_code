function [cohesiveZone, T_mean,PhediLoc_mean, newPowerStruct] = phedi_CohesiveFromFit(DataStruct,varargin)
% cohesiveZone = phedi_CohesiveFromFit(DataStruct, Cf, t_min, t_max, Shifts, selectedPhedis, doPlot)
%
% only DataStruct is mandatory
%% parameters
%--- set defaults
Cf_def = DataStruct.PhediData.Cf;
[Cf, t_min, t_max, Shifts, selectedPhedis, doPlot] = setDefaults4function(varargin,...
    Cf_def,...
    0.003/Cf_def,...
    0.025/Cf_def,...
    linspace(-1.5*1e-4,1.5*1e-4,200),...
    find(DataStruct.PhediData.slopeIncline==1),...
    0);
if isempty(varargin{2})
    t_min = 0.003/Cf;
end
if isempty(varargin{3})
    t_max = 0.025/Cf;
end

%% analyze power
[T_mean,PhediLoc_mean] = phedi_averagePhedi(DataStruct, selectedPhedis,[],[],0);
fitPowerStruct = phedi_fitPowerToLocation(T_mean,PhediLoc_mean,t_min, t_max,Shifts);

%% find 0.5 fit
[~,divergeIdx] = min(abs(fitPowerStruct.MDLs(1,:)-0.5));
[relevantShift,halfPower] = intersections(fitPowerStruct.Shifts,fitPowerStruct.MDLs(1,:),[-1 1]*1e9,0.5*[1 1]);
relevantShift = mean(relevantShift);
newPowerStruct = phedi_fitPowerToLocation(T_mean,PhediLoc_mean,t_min, t_max,relevantShift);
try
XX = T_mean+relevantShift;
catch
    disp();
end
YY_fit = (10^newPowerStruct.MDLs(2))*(XX.^0.5);
YY_fit(XX<=0) = 0;

%% find cohesive
thresh = 1e-7;
modelDiff = PhediLoc_mean - YY_fit;
inModelLogical = abs(modelDiff)<thresh;

%--- find zero nearest to origin
[~,nearZeroIdx] = min(abs(XX(~inModelLogical)));
CC = bwlabel(~inModelLogical);
selectedCC = CC(~inModelLogical);
cohesiveZoneVec = XX(CC==selectedCC(nearZeroIdx));
cohesiveZone = [min(cohesiveZoneVec),max(cohesiveZoneVec)];

%% plot stuff
if doPlot
    figure;
    plot(-T_mean*Cf,PhediLoc_mean);
    hold on;
    cohesive_logic = T_mean>min(cohesiveZone) & T_mean<max(cohesiveZone);
    plot(-Cf*T_mean(cohesive_logic),PhediLoc_mean(cohesive_logic),'r','LineWidth',1.2);
    plot(-Cf*XX,YY_fit,'k','LineWidth',1.2);
    xlabel('x-x_{tip} [m]');
    ylabel('displacement [m]');
    
end

end