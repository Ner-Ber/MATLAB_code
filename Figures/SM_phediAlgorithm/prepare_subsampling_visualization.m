%% Make a simultion of phedis and tracking
doSimulate = 0;
if doSimulate
    % close all;
    % clear;
    
    %% set "experimetn" parameters
    %-- camera specs and setup settings
    fps = 1e6;
    res = 97500;
    FillFac = 56/100; % in fraction of surface
    scratchWidth = 0.0001;
    simpleSamp = 0;
    frameLength = 384;
    
    %-- algorithm specs
    prePhediTime = 0.3e-4;
    postPhediTime = 5e-4;
    numOOM = 3;
    smooth4Phedi = 0.001;
    PhediSafetyDist = 160;
    
    relative2first = 0;
    coarseShift = 0;
    % if relative2first % if relative 2 first impose coarse shift
    %     coarseShift = 1;
    % end
    QuickAndDirty = 1;
    P = [0.1 0.375 0.5 0.75 0.9 1];
    interpP = P(3);
    slopePenalty = -0;
    expandLowBound = -1;
    expandHighBound = 0;
    
    Cf = 120;
    %% create RowOverTime
    %--- generate "true" profile
    % dx = 0.01;
    dx = 0.05;
    x_true = 0.5:dx:(frameLength+0.5);                % sub pixles coordinates
    % x_true = 1:dx:frameLength;                % sub pixles coordinates
    res_true = res/dx;
    spatialVec_true = x_true./res;
    
    %-- sin:
    if 0
        intns_vec_true = sin(pi*spatialVec_true/(2*scratchWidth)).^2;
        %     intns_vec_true = sin(x_true*2*pi/19).^2;
    end
    
    %-- assymetric betafunc:
    if 0
        y1 = betapdf(mod(x_true/12,1.5),5,2);
        intns_vec_true = y1/max(y1);
        intns_vec_true = (intns_vec_true);  % modification for checking
    end
    
    %-- Gaussians
    if 0
        stepWidth = 1.5e-1;
        G = @(X) exp(-((mod(X/2/scratchWidth,1)-0.5)/stepWidth).^2);
        intns_vec_true = G(spatialVec_true);
    end
    
    %--- smoothed step
    if 0
        S_W = [0.005 0.025 0.05 0.125 0.18];
        stepWidth = S_W(1);
        S = @(x) ((tanh((mod(x/2/scratchWidth,1)-0.25)./stepWidth).*tanh((0.75-mod(x/2/scratchWidth,1))./stepWidth))+1).*0.5;
        intns_vec_true = S(spatialVec_true);
    end
    
    %--- interference pattern
    if 1
        Y = 0.002;
        stepWidth = 0.003;
        Lambda = 532*1e-9;
        [intns_vec_true, Psi0] = FresnelCreateSlitsPattern(spatialVec_true, Y, scratchWidth, stepWidth, Lambda);
        intns_vec_true = abs(intns_vec_true).^2;
    end
    
    %--- create time vector
    Nframes = 6e3;
    t_vec = (-Nframes/2:Nframes/2);
    t_vec_s = t_vec/fps;
    zeroFrame = find(t_vec==0);
    %--- create an unshifted RowOverTime
    RowOverTime_true_unshifted = repmat(intns_vec_true,length(t_vec),1);
    
    %% shift row over time
    %--- simple shift
    if 1
        if 0
            XX = sort(-t_vec_s*Cf);
            [Cd, Cs, Cr, nu, ro, E, mu, Gamma, PlaneStrain, tau_p, Xc0]=CrackSolutionMaterialProperties;
            warning('off');
            % Sol_cohsvMdl = CrackSolutionGeneralCohesive_InsertVariables_Neri(Cf/Cr,2e-9,6,'L',2.5e-3,'const',[],[min(XX) max(XX)],XX(2)-XX(1));
            Sol_cohsvMdl = CrackSolutionGeneralCohesive_InsertVariables_Neri(Cf/Cr,2e-9,6,'L',2.5e-3,'const',[],rot90(XX,2));
            warning('on');
            shiftVector = round(Sol_cohsvMdl.Ux*res_true);
        end
        
        %-- constant slifind profile
        if 1
            Vslip = Cf/700;
            Vslip_true_pix_per_frame = (Vslip*res_true)./fps;
            shiftVector = zeros(size(t_vec));
            shiftVector(zeroFrame:end) = (1:(length(shiftVector)-(zeroFrame-1)));
            shiftVector = round(Vslip_true_pix_per_frame*shiftVector);
            shiftVector(zeroFrame:end) = Vslip_true_pix_per_frame;
            shiftVector = round(cumsum(shiftVector));
        end
        
        
        % shiftVector = 0.5*t_vec+5*sin(0.1*t_vec);
        % shiftVector(t_vec<=0) = 0;
        % shiftVector(t_vec>=300) = shiftVector(t_vec==300);
        % shiftVector = round(shiftVector);
        % shiftVector = -shiftVector;         % modification for checking
        %--- create shifted RowOverTime
        RowOverTime_true = MyCircshift(RowOverTime_true_unshifted,shiftVector);
    end
    
    %--- shift with strains
    if 0
        Cf_PIXperFram = Cf.*res_true./fps;
        slope = -1/Cf_PIXperFram ;
        
        
    end
    
    %% create a measurement row over time
    %--- create valiable RowOverTime
    [RowOverTime_preDamp, sampleLocations] = phediSimulation_sampleFromRealROT(...
        RowOverTime_true, x_true, frameLength, simpleSamp, FillFac);
    % sampleLocations = ismember(x_true,1:frameLength); %~mod(x_true,1); % find pixle locations
    % if simpleSamp   % this case is just undersampling the data on the x direction
    %     RowOverTime_preDamp = RowOverTime_true(:,sampleLocations);
    % else    % this option will integrate by a given fill factor
    %     x_true_mod = mod(x_true,1);
    %     D = abs(sqrt(FillFac))/2;
    %     fillLogic = x_true_mod<=D | x_true_mod>=(1-D);
    %     RowOverTime_fill = RowOverTime_true(:,fillLogic);
    %     ROT_reshape = reshape(RowOverTime_fill,size(RowOverTime_fill,1),[],frameLength);
    %     RowOverTime_preDamp = permute(sum(ROT_reshape,2),[1,3,2])./(size(ROT_reshape,2));
    % end
    
    %--- create amplitude decrease
    RowOverTime = RowOverTime_preDamp;
    %-- step decrease
    if 1
        RowOverTime(t_vec>0,:) = 0.8*RowOverTime(t_vec>0,:);
    end
    %-- exponential decrease
    if 1
        tau = 10; % in frames
        decay = exp(-t_vec./tau)*0.2 + 0.8;
        decay(t_vec<=0) = 1;
        decay = decay(:);
        RowOverTime = bsxfun(@times, RowOverTime, decay);
    end
    
    %--- create stronger light in center
    x0 = mean(x_true);
    GaussFilt = @(x) exp(-((x-x0)./(length(x_true)/85)).^2);
    GaussMask = GaussFilt(x_true(sampleLocations))./max(x_true(sampleLocations));
    RowOverTime = bsxfun(@times,RowOverTime,GaussMask);
    
    %--- blur RowOverTime
    BlurAmount = 9;
    % g = createGaussianAproxCostume(1,BlurAmount);
    g = createGaussianAproxCostume(BlurAmount);
    % RowOverTime = conv2(RowOverTime,g,'same');
    
    %--- add noise ROT
    Noise = randn(size(RowOverTime));
    Noise = Noise./max(Noise(:));
    RowOverTime = RowOverTime+Noise*0; %1e-4;
    
    
    %% follow and analyze asperities
    disp('Definning and following asperities...');
    [vallyPairs,~,maxTrajectories,CM_Trajectories,skewnessMat,varMat,center_of_mass, MarkedMatCell]...
        = Movie_followAsperities(RowOverTime);
    spatialVec = (0:(size(RowOverTime,2)-1))/res;
    
    
    %% calculate phedis
    disp('analyzing phedis...');
    RowOverTime_normalized = ROT_normalizeIntesByAsperity(RowOverTime,MarkedMatCell);
    
    %--- define frame to calc phedis on:
    prePhediFrames = round(prePhediTime*fps);
    postPhediFrames = round(postPhediTime*fps);
    firstFrame4Phedi = Nframes/2-prePhediFrames;
    lastFrame4Phedi = Nframes/2+postPhediFrames;
    if firstFrame4Phedi<1; firstFrame4Phedi=1; end;
    if lastFrame4Phedi>size(RowOverTime_normalized,1); lastFrame4Phedi=size(RowOverTime_normalized,1); end;
    
    %--- get phedis initial locations
    FirstRowSig = RowOverTime_normalized(firstFrame4Phedi,:);
    phedisInitialLocsInPix = phedi_getPhediInitialLocs(FirstRowSig,smooth4Phedi, PhediSafetyDist, res, scratchWidth);
    phedisInitialLocsInPix = phedisInitialLocsInPix(1:2);
    
    
    %--- calculate coarse shift
    if coarseShift
        [coarseShiftsVector, shiftedRowOverTime] = Movie_findCoarseShiftsInRowOverTime(RowOverTime_normalized,center_of_mass);
    else
        shiftedRowOverTime = RowOverTime_normalized;
        coarseShiftsVector = zeros(size(center_of_mass,1),1);
    end
    %--- follow phedis
    calShiftTimer = tic;
    %--- total function
    try
        [frameCount, PhediLocationPixShifted, slopeInclineMat, StrechingFactorsMat, SlopeRegions] = Movie_phedi_calculate_continues_shift_total(...
            shiftedRowOverTime, firstFrame4Phedi, lastFrame4Phedi, phedisInitialLocsInPix,...
            'RMS_minimization',numOOM, slopePenalty, expandLowBound, expandHighBound, relative2first);
    catch
        disp('');
    end
    
    
    %--- fit with previous single frame
    % [frameCount, PhediLocationPixShifted, slopeIncline, StrechingFactorsMat, SlopeRegions] = Movie_phedi_calculate_continues_shift(...
    %     shiftedRowOverTime, firstFrame4Phedi, lastFrame4Phedi, phedisInitialLocsInPix,numOOM,phediRelative2first,slopePenalty, QuickAndDirty,[],interpP);
    
    %--- fit with fresnel function
    % [frameCount, PhediLocationPixShifted, slopeIncline, StrechingFactorsMat, strechFactorHorzMat, ShiftsMat] =...
    %     Movie_phedi_calcShiftByFit(...
    %     shiftedRowOverTime, firstFrame4Phedi, lastFrame4Phedi, phedisInitialLocsInPix);
    
    %--- fit with fresnel function using previous frames
    % [frameCount, PhediLocationPixShifted, slopeIncline, VertStrechFactorMat, HorzSrechFactorMat, HorzDisplacementMat] =...
    %         Movie_phedi_calcShiftByFit_wPerturbtion(...
    %         shiftedRowOverTime, firstFrame4Phedi, lastFrame4Phedi, phedisInitialLocsInPix);
    
    relevantCoarseShifts = coarseShiftsVector(firstFrame4Phedi:lastFrame4Phedi);
    if nnz(coarseShiftsVector>0)
        PhediLocationPix = bsxfun(@plus,PhediLocationPixShifted,relevantCoarseShifts)-min(coarseShiftsVector(coarseShiftsVector>0));
    else
        PhediLocationPix = bsxfun(@plus,PhediLocationPixShifted,relevantCoarseShifts);
    end
    
    
    %--- transform units: (minus 1 since the meters scale starts from 0 and the
    %pixel scale starts from 1)
    PhediLocation = (PhediLocationPix-1)/res;
    measuredPhedisFromPlot = spatialVec(phedisInitialLocsInPix);
    timeVec= t_vec_s(frameCount);
    
    
    T_shiftTimer = toc(calShiftTimer);
    disp(['time for calculating phedi shifts = ',num2str(T_shiftTimer)]);
    
    %%
    slopeIncline = slopeInclineMat(1,:);
    % slopeIncline = mean(slopeInclineMat);
    
    
    slopeIncline_originalAll = slopeIncline;
    slopeIncline = slopeIncline(1,:);
    
    %% plot results
    slope2Plot = 0;
    plotDamaged = 1;
    Nphd = size(PhediLocation,2);       % total number of phedis
    FigColors = MyVaryColor(Nphd);
    FigColors(slopeIncline==0,:) = repmat([1 0 0],nnz(slopeIncline==0),1);    % set damages to red
    %--names:
    LGND = cellfun(@num2str,num2cell(measuredPhedisFromPlot),'UniformOutput',0);
    %--markers:
    allMarks = {'*-','.-','--'};
    marksLookup = slopeIncline;
    marksLookup(slopeIncline==-1) = 2;
    marksLookup(slopeIncline==0) = 3;
    FigMarks = allMarks(marksLookup);
    %--phedis to plot by number, pick relevant phedis index
    phdIndx = [];
    if plotDamaged
        phdIndx = cat(2,phdIndx,find(slopeIncline==0));
    end
    if slope2Plot==-1
        phdIndx = cat(2,phdIndx,find(slopeIncline==-1));
    elseif slope2Plot==1
        phdIndx = cat(2,phdIndx,find(slopeIncline==1));
    else
        phdIndx = cat(2,phdIndx,find(slopeIncline~=0));
    end
    phdIndx = unique(phdIndx);
    
    
end

%% make a figure of phedis
ROT_chosen = RowOverTime_preDamp;
figure;
% subplot(1,2,1);
plot(spatialVec, ROT_chosen(firstFrame4Phedi,:),'r*','MarkerSize',6);
hold on;
plot((x_true-1)/res, RowOverTime_true(firstFrame4Phedi,:),'b-');
% plot(spatialVec(phedisInitialLocsInPix),RowOverTime_preDamp(firstFrame4Phedi,phedisInitialLocsInPix),'ro');

N = 80;
plot((x_true-1)/res, RowOverTime_true(firstFrame4Phedi+N,:),'.-','color',rgb('DarkCyan'),'DisplayName',num2str(N));
plot(spatialVec, ROT_chosen(firstFrame4Phedi+N,:),'*','MarkerSize',6, 'DisplayName',num2str(N),'color',rgb('IndianRed'));
N=120;
plot((x_true-1)/res, RowOverTime_true(firstFrame4Phedi+N,:),'.-','color',rgb('PowderBlue'),'DisplayName',num2str(N));
plot(spatialVec, ROT_chosen(firstFrame4Phedi+N,:),'*','MarkerSize',6, 'DisplayName',num2str(N),'color',rgb('LightSalmon'));

xlabel('aprox. location on block');
ylabel('normalized  intensity [arb]')
