%% follow phedi simulation
% close all;
% clear;

fps = 1e6;
prePhediTime = 1.5e-4;
postPhediTime = 1.5e-4;
numOOM = 2;
%% create RowOverTime
%--- generate "true" profile
dx = 0.25;
x_true = 1:dx:384;                % sub pixles coordinates
% intns_vec_true = sin(x_true*2*pi/19).^2;
y1 = betapdf(mod(x_true/12,1.5),5,2);
intns_vec_true = y1/max(y1);
intns_vec_true = fliplr(intns_vec_true);  % modification for checking

%--- create time vector
Nframes = 6e3;
t_vec = -Nframes/2:Nframes/2;
zeroFrame = find(t_vec==0);
%--- create an unshifted RowOverTime
RowOverTime_true_unshifted = repmat(intns_vec_true,length(t_vec),1);
%--- create the shifting (sliding) profile
shiftVector = 0.5*t_vec+5*sin(0.1*t_vec);
shiftVector(t_vec<=0) = 0;
shiftVector(t_vec>=300) = shiftVector(t_vec==300);
shiftVector = round(shiftVector);
shiftVector = -shiftVector;         % modification for checking
%--- create shifted RowOverTime
RowOverTime_true = MyCircshift(RowOverTime_true_unshifted,shiftVector);
%--- create valiableRowOverTime
sampleLocations = ~mod(x_true,1); % find pixle locations
RowOverTime = RowOverTime_true(:,sampleLocations);
%--- create amplitude decrease
RowOverTime(t_vec>0,:) = 0.8*RowOverTime(t_vec>0,:);
%--- create stronger light in center
x0 = mean(x_true);
GaussFilt = @(x) exp(-((x-x0)./(length(x_true)/20)).^2);
GaussMask = GaussFilt(x_true(sampleLocations))./max(x_true(sampleLocations));
RowOverTime = bsxfun(@times,RowOverTime,GaussMask);
%--- add noise ROT
Noise = randn(size(RowOverTime));
Noise = Noise./max(Noise(:));
RowOverTime = RowOverTime+Noise*1e-4;

%% find center of mass trajectory
[vallyPairsCell,~] = Movie_define_asperities_for_scratches({RowOverTime});
[~,~, ~,MarkedMatCell, ~] = Movie_follow_maxima_in_row(permute(RowOverTime,[3 2 1]),...
    vallyPairsCell, 1, 'all', [], 1, [],'min');
[~,CM_Trajectories,~,~,center_of_massCell] = Movie_analyze_asperities(MarkedMatCell, {RowOverTime});

%% calculate phedis
disp('analyzing phedis...');
relevanrRowOverTime = RowOverTime;

%--- define frame to calc phedis on:
prePhediFrames = round(prePhediTime*fps);
postPhediFrames = round(postPhediTime*fps);
firstFrame4Phedi = Nframes/2-prePhediFrames;
lastFrame4Phedi = Nframes/2+postPhediFrames;

smoothParameter = 0.001;
phedisInitialLocsInPix = Movie_phedi_findTrenches(relevanrRowOverTime(firstFrame4Phedi,:), smoothParameter);
%--- eliminate phedis close to edge:
safetyDist = 70;
phedisInitialLocsInPix = phedisInitialLocsInPix(abs(phedisInitialLocsInPix-length(relevanrRowOverTime(firstFrame4Phedi,:)))>safetyDist);
phedisInitialLocsInPix = phedisInitialLocsInPix(phedisInitialLocsInPix>safetyDist);
phedisInitialLocsInPix = phedisInitialLocsInPix(2:(end-1));

useShiftedMethod = 1;
if useShiftedMethod
    %--- calculate coarse shift
    [coarseShiftsVector, shiftedRowOverTime] = Movie_findCoarseShiftsInRowOverTime(relevanrRowOverTime,center_of_massCell{end});
    relevanrRowOverTime = shiftedRowOverTime;
end


calShiftTimer = tic;
if useShiftedMethod
    [frameCount, PhediLocationPixShifted, slopeIncline, StrechingFactorsMat] = Movie_phedi_calculate_continues_shift(...
        relevanrRowOverTime, firstFrame4Phedi, lastFrame4Phedi, phedisInitialLocsInPix,numOOM,1);
    relevantCoarseShifts = coarseShiftsVector(firstFrame4Phedi:lastFrame4Phedi);
    try
        PhediLocationPix = bsxfun(@plus,PhediLocationPixShifted,relevantCoarseShifts)-min(coarseShiftsVector(coarseShiftsVector>0));
    catch
        PhediLocationPix = bsxfun(@plus,PhediLocationPixShifted,relevantCoarseShifts)-max(coarseShiftsVector(coarseShiftsVector<0));
    end
    
    
else
    %     [frameCount, PhediLocationPix, slopeIncline, StrechingFactorsMat] = Movie_phedi_calculate_continues_shift(...
    %         relevanrRowOverTime, firstFrame4Phedi, lastFrame4Phedi, phedisInitialLocsInPix,numOOM,0);
    %     [frameCount, PhediLocationPix, slopeIncline, StrechingFactorsMat] = Movie_phedi_calculate_continues_shift_updating(...
    %         relevanrRowOverTime, firstFrame4Phedi, lastFrame4Phedi, phedisInitialLocsInPix,numOOM,0);
end
T_shiftTimer = toc(calShiftTimer);
disp(['time for calculating phedi shifts = ',num2str(T_shiftTimer)]);

%% plot results
%--- plot first frame for phedi
figure;
plot(RowOverTime(firstFrame4Phedi,:));
hold on;
plot(phedisInitialLocsInPix,RowOverTime(firstFrame4Phedi,phedisInitialLocsInPix),'ro');
legend('first frame of phedi analysis','phedi locations');
hold off;
title('first row of phedi analysis');

%--- plot rowOverTime and phedis
figure;
PhediTrajectories = bsxfun(@plus, PhediLocationPix,phedisInitialLocsInPix(:)');
imagesc(RowOverTime);
hold on;
plot(PhediTrajectories, firstFrame4Phedi:lastFrame4Phedi,'LineWidth',1.5);
hold off;
title('ROT with phedis');

%--- plot original translation and phedi location
figure;
plot(t_vec(firstFrame4Phedi:lastFrame4Phedi),shiftVector(firstFrame4Phedi:lastFrame4Phedi)*dx,'.-');
hold on;
plot((firstFrame4Phedi:lastFrame4Phedi)-zeroFrame,...
    bsxfun(@minus,PhediLocationPix(:,slopeIncline==1),PhediLocationPix(1,slopeIncline==1)),...
    '*');
plot((firstFrame4Phedi:lastFrame4Phedi)-zeroFrame,...
    bsxfun(@minus,PhediLocationPix(:,slopeIncline==-1),PhediLocationPix(1,slopeIncline==-1)),...
    'd');
title({'original translation and phedi trajectory','.=original   *=left side   d=right side'});