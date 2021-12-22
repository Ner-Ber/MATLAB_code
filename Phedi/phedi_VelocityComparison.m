function [rnkMaxCell, rnkP2pCell, maxVelCell, p2pVelCell] = phedi_VelocityComparison(PhediStructCell,varargin)
%[rnkMaxCell, rnkP2pCell, maxVelCell, p2pVelCell] = phedi_VelocityComparison(PhediStructCell,DoPlot)
%
%  OUTPUTS:
%rnkMaxCell - ranking of velocities amplitudes from highest to lowest in each event.
%rnkP2pCell - ranking of velocities oscillation amps from highest to lowest in each event.
%maxVelCell - values of the velocity peaks in each phedi in each event
%p2pVelCell - values of the velocity oscillation amp in each phedi in each event
%
% INPUTS:
%DoPlot (optional): default 0


[DoPlot] = setDefaults4function(varargin,0);

relevantEvents = find(~cellfun(@isempty,PhediStructCell) &...
    ~cellfun(@(a) isfield(a,'theWorkspace'),PhediStructCell));

if DoPlot
    N = length(relevantEvents);
    n = ceil(3*N/12);
    m = ceil(N/n);
    figure;
    ha = tight_subplot(n,m,[.02 .03],[.1 .01],[.01 .01]);
end

time_interval = [-0.5 0.5]*1e-4;
maxVelCell = {};
p2pVelCell = {};
rnkMaxCell = {};
rnkP2pCell = {};
for i=relevantEvents(:)'
    %--- cehck whether the event happened (exclude precursurs)
    scratchRegionMeters = PhediStructCell{i}.ExperimentData.scratchRegionMeters;
    frontFoundX = PhediStructCell{i}.BigPicRotStruct.x(PhediStructCell{i}.BigPicRotStruct.frontStepsPix);
    stepsInScratches = nnz(frontFoundX<max(scratchRegionMeters) & frontFoundX>min(scratchRegionMeters));
    if stepsInScratches<2
        continue
    end
    
    chckVelLogical = ...
        PhediStructCell{i}.PhediData.t_mins_t_tip_4vel>min(time_interval) &...
        PhediStructCell{i}.PhediData.t_mins_t_tip_4vel<max(time_interval);
    %-- iterate over phedis:
    maxVelVec = [];
    p2pVelVec = [];
    maxVelTime = [];
    for phediIdx = find(~~sum(chckVelLogical))
        relevantVel = PhediStructCell{i}.PhediData.PhediVelocity(chckVelLogical(:,phediIdx),phediIdx);
        relevantTime = PhediStructCell{i}.PhediData.t_mins_t_tip_4vel(chckVelLogical(:,phediIdx),phediIdx);
        [thisMaxVel,VelIdx] = max(relevantVel);
        maxVelTime(phediIdx) = relevantTime(VelIdx);
        maxVelVec(phediIdx)= thisMaxVel;
        p2pVelVec(phediIdx) = peak2peak(relevantVel);
    end
    if DoPlot
        axes(ha(i==relevantEvents));
%         figure;
        hold on;
        plot(PhediStructCell{i}.PhediData.t_mins_t_tip_4vel,PhediStructCell{i}.PhediData.PhediVelocity,'.-');
        plot(maxVelTime,maxVelVec,'ro');
        title(['Ev=',num2str(i)]);
    end
    [~,~,rnkMax] = unique(maxVelVec);
    [~,~,rnkP2p] = unique(p2pVelVec);
    rnkMaxCell{i} = rnkMax;
    rnkP2pCell{i} = rnkP2p;
    maxVelCell{i} = maxVelVec;
    p2pVelCell{i} = p2pVelVec;
end


end