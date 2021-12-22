%% define parameters
Lambda = 532*1e-9;
k = 2*pi/Lambda;
res = 97500;
dx = 0.08;
frameLength = 384;
x_true = 0.5:dx:(frameLength+0.5);                % sub pixles coordinates
res_true = res/dx;

scratchWidth = 0.0001;
S_W = [0.005 0.025 0.05 0.125 0.18];
stepWidth = S_W(1);
%% prepare initial conditions
Nsteps = 11;
% distReativeD = 1:-0.02:0.5;
distReativeD = 0:2e-6:2e-5;

spatialVec_true = -0.6e-3:8e-7:0.6e-3;
I = nan(length(distReativeD),length(spatialVec_true));
PsiMat = nan(length(distReativeD),length(spatialVec_true));
Ishifted = I;
Psi0 = @(p) ((tanh((p/2/scratchWidth+0.25)./stepWidth).*tanh((0.25-p/2/scratchWidth)./stepWidth))+1).*0.5;

for j = 1:length(distReativeD)
%     D = (Nsteps-1)*scratchWidth*distReativeD(j);
    D = (Nsteps-1)*scratchWidth;
    Dunit = 2*D/(Nsteps-1);
    
    combLocation = -D:Dunit:D;
    combLocation((Nsteps+1)/2+1) = combLocation((Nsteps+1)/2+1) - distReativeD(j);
    
    [~,combLocationPix] = min(abs(bsxfun(@minus,spatialVec_true(:),combLocation(:)')),[],1);
    comb = zeros(size(spatialVec_true));
    comb(combLocationPix) = 1;
    Psi0conv = @(p) conv(Psi0(p),comb,'same');
    Psi0Damped = @(p) Psi0conv(p);
    
    %% solve integral
    Y = 0.002;
    warning('off');
    yj = Y;
    for i = 1:length(spatialVec_true);
        x_i = spatialVec_true(i);
        Intgrnd = @(p) Psi0Damped(spatialVec_true).*sqrt(k./(2*pi*1i*yj)).*exp(1i.*k.*(yj+(x_i-spatialVec_true).^2./(2*yj)));
        I(j,i) = trapz(spatialVec_true,Intgrnd(spatialVec_true));
    end
    Ishifted(j,:) = I(j,:);
    PsiMat(j,:) = Psi0Damped(spatialVec_true);
    warning('on');
end

%% single step
PsiStep = @(p) ((tanh((0.25-p/2/scratchWidth)./stepWidth))+1).*0.5;
Istep = nan(1,length(spatialVec_true));
for i = 1:length(spatialVec_true);
    x_i = spatialVec_true(i);
    Intgrnd = @(p) PsiStep(spatialVec_true).*sqrt(k./(2*pi*1i*yj)).*exp(1i.*k.*(yj+(x_i-spatialVec_true).^2./(2*yj)));
    Istep(i) = trapz(spatialVec_true,Intgrnd(spatialVec_true));
end
Ishifted = [Istep;Ishifted];

%% create difference mat
IshiftedDiffMat = diff(abs(Ishifted),1,1);
relevantDiffRegionLogic = spatialVec_true>=-scratchWidth & spatialVec_true<=scratchWidth;

IshiftedRatioMat = abs(Ishifted(2:end,:))./abs(Ishifted(1:end-1,:));


%% plot
figure;
hold on;
plot(spatialVec_true,abs(Ishifted'),'LineWidth',1.5);
plot(spatialVec_true,Psi0(spatialVec_true),'DisplayName','original');
plot(spatialVec_true,PsiStep(spatialVec_true),'DisplayName','original step');
title('interference patterns');

figure; hold on;
subplot(2,1,1); hold on;
plot(spatialVec_true(relevantDiffRegionLogic),IshiftedDiffMat(:,relevantDiffRegionLogic),'LineWidth',1.5);
xlim([0 1]*1e-4);
legend({'step-1','2-1','3-2','4-3','5-4','6-5','7-6'});
title('difference between following scenarios');
subplot(2,1,2); hold on;
plot(spatialVec_true(relevantDiffRegionLogic),IshiftedDiffMat(:,relevantDiffRegionLogic),'LineWidth',1.5);
xlim([0 1]*1e-4);


figure; hold on;
plot(spatialVec_true(relevantDiffRegionLogic),IshiftedRatioMat(:,relevantDiffRegionLogic),'LineWidth',1.5);
xlim([0 1]*1e-4);
legend({'step-1','2-1','3-2','4-3','5-4','6-5','7-6'});
title('ratio between following scenarios');