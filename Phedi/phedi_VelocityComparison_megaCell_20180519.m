function [rnkMaxCell, rnkP2pCell, maxVelCell, p2pVelCell, stretchCell, rnkStrchCell, shiftsCell, rnkShiftCell,...
        GammaCell, LocalGammaCell, rnkLocalGCell]...
        = phedi_VelocityComparison_megaCell_20180519(PhediStructCellOfCell,varargin)
    %[rnkMaxCell, rnkP2pCell, maxVelCell, p2pVelCell, stretchCell, GammaCell] = phedi_VelocityComparison(PhediStructCell,DoPlot)
    %
    %  OUTPUTS:
    %rnkMaxCell - ranking of velocities amplitudes from highest to lowest in each event.
    %rnkP2pCell - ranking of velocities oscillation amps from highest to lowest in each event.
    %maxVelCell - values of the velocity peaks in each phedi in each event
    %p2pVelCell - values of the velocity oscillation amp in each phedi in each event
    %
    % INPUTS:
    %DoPlot (optional): default 0
    %
    %
    % THIS IS AN OLD VERSION THAT TAKES A 'CELL OF CELLS' WITH ALL EVENTS
    % (RELEVANT AND NON) AND FILTERS THEM HERE INSIDE THE FUNCTION.
    % THE NEWER VERSION WILL TAKE A FILTERED 'CELL' (ONLY).
    %
    
    
    [DoPlot] = setDefaults4function(varargin,0);
    
    PhediStructCellMega = [PhediStructCellOfCell{:}];
    % relevantEvents = find(~cellfun(@isempty,PhediStructCellMega) & ~cellfun(@(a) isfield(a,'theWorkspace'),PhediStructCellMega));
    relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCellMega));
    
    % relevantEvents = find(~cellfun(@isempty,PhediStructCell) &...
    %     ~cellfun(@(a) isfield(a,'theWorkspace'),PhediStructCell));
    [Cd, Cs, Cr, nu, ro, E, mu, Gamma, PlaneStrain, tau_p, Xc0]=CrackSolutionMaterialProperties;
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
    rnkStrchCell = {};
    stretchCell = {};
    shiftsCell = {};
    rnkShiftCell = {};
    GammaCell = {};
    LocalGammaCell = {};
    rnkLocalGCell = {};
    for i=relevantEvents(:)'
        %--- cehck whether the event happened (exclude precursurs)
        scratchRegionMeters = PhediStructCellMega{i}.ExperimentData.scratchRegionMeters;
        frontFoundX = PhediStructCellMega{i}.BigPicRotStruct.x(PhediStructCellMega{i}.BigPicRotStruct.frontStepsPix);
        stepsInScratches = nnz(frontFoundX<max(scratchRegionMeters) & frontFoundX>min(scratchRegionMeters));
        if stepsInScratches<2
            continue
        end
        
        chckVelLogical = ...
            PhediStructCellMega{i}.PhediData.t_mins_t_tip_4vel>min(time_interval) &...
            PhediStructCellMega{i}.PhediData.t_mins_t_tip_4vel<max(time_interval);
        [~,~, stretch_vec] = phedi_UxBestFitVerticalStretch(PhediStructCellMega{i});
        [~,~, shifts_vec] = phedi_UxBestFitVerticalShift(PhediStructCellMega{i});
        %-- iterate over phedis:
        maxVelVec = [];
        p2pVelVec = [];
        maxVelTime = [];
        stretchVec = [];
        shiftsVec = [];
        GammaVec = [];
        LocalGammaVec = [];
        for phediIdx = find(~~sum(chckVelLogical))
            relevantVel = PhediStructCellMega{i}.PhediData.PhediVelocity(chckVelLogical(:,phediIdx),phediIdx);
            relevantTime = PhediStructCellMega{i}.PhediData.t_mins_t_tip_4vel(chckVelLogical(:,phediIdx),phediIdx);
            [thisMaxVel,VelIdx] = max(relevantVel);
            maxVelTime(phediIdx) = relevantTime(VelIdx);
            maxVelVec(phediIdx)= thisMaxVel;
            p2pVelVec(phediIdx) = peak2peak(relevantVel);
            stretchVec(phediIdx) = stretch_vec(phediIdx);
            shiftsVec(phediIdx) = shifts_vec(phediIdx);
            GammaVec(phediIdx) = PhediStructCellMega{i}.solAtSG.Gamma;
            
            %-- local gamma vec
            fixedK = PhediStructCellMega{i}.solAtSG.K*stretch_vec(phediIdx);
            v = PhediStructCellMega{i}.PhediData.Cf;
            k=Cs/Cd;%Broberg p.330
            alpha_d=(1-(v./Cd).^2).^0.5;
            alpha_s=(1-(v./Cs).^2).^0.5;
            D=4*alpha_d.*alpha_s-(1+alpha_s.^2).^2;
            %----- Following Broberg p.334,336
            A=2*(1-k^2)*alpha_s.*v.^2./D/Cs^2;
            localG = (fixedK^2)./(mu*4*(1-k^2)./A);
            LocalGammaVec(phediIdx) = localG;
        end
        if DoPlot
            axes(ha(i==relevantEvents));
            %         figure;
            hold on;
            plot(PhediStructCellMega{i}.PhediData.t_mins_t_tip_4vel,PhediStructCellMega{i}.PhediData.PhediVelocity,'.-');
            plot(maxVelTime,maxVelVec,'ro');
            title(['Ev=',num2str(i)]);
        end
        [~,~,rnkMax] = unique(maxVelVec);
        [~,~,rnkP2p] = unique(p2pVelVec);
        [~,~,rnkStrtch] = unique(stretchVec);
        [~,~,rnkShift] = unique(shiftsVec);
        [~,~,rnkLoclG] = unique(LocalGammaVec);
        rnkMaxCell{i} = rnkMax;
        rnkP2pCell{i} = rnkP2p;
        maxVelCell{i} = maxVelVec;
        p2pVelCell{i} = p2pVelVec;
        stretchCell{i} = stretchVec;
        shiftsCell{i} = shiftsVec;
        rnkStrchCell{i} = rnkStrtch;
        rnkShiftCell{i} = rnkShift;
        GammaCell{i} = GammaVec;
        LocalGammaCell{i} = LocalGammaVec;
        rnkLocalGCell{i} = rnkLoclG;
    end
    
    
end