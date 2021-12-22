%% batch run phedi simulation
% close all

if 0
    P = [0.1 0.375 0.5 0.75 0.9 1];
    S_W = [0.005 0.025 0.05 0.125 0.18];
    numOOM_vec = 2:3;
    
    phediDiff_cell = cell(length(P),length(S_W),length(numOOM_vec));
    
    driftVel1_MAT = nan(length(P),length(S_W),length(numOOM_vec));
    driftVel2_MAT = nan(length(P),length(S_W),length(numOOM_vec));
    
    oscilT1_MAT = nan(length(P),length(S_W),length(numOOM_vec));
    oscilT2_MAT = nan(length(P),length(S_W),length(numOOM_vec));
    
    Amp1_MAT = nan(length(P),length(S_W),length(numOOM_vec));
    Amp2_MAT = nan(length(P),length(S_W),length(numOOM_vec));
    
    for i_p = 1:length(P)
        for i_sw = 1:length(S_W)
            for i_oom = 1:length(numOOM_vec)
                interpP = P(i_p)
                stepWidth = S_W(i_sw)
                numOOM = numOOM_vec(i_oom)
                
                follow_phedi_simulation;
                phediDiff_cell{i_p,i_sw,i_oom} = phediDiff;
                phediDiff_smth1 = smooth(phediDiff(timeVec>0,1),11,'lowess');
                phediDiff_smth2 = smooth(phediDiff(timeVec>0,2),11,'lowess');
                [pks1,locs1] = findpeaks(phediDiff_smth1,timeVec(timeVec>0));
                [pks2,locs2] = findpeaks(phediDiff_smth2,timeVec(timeVec>0));
                %-- drift vel
                driftVel1_MAT(i_p,i_sw,i_oom) = mean(diff(pks1(:))./diff(locs1(:)));
                driftVel2_MAT(i_p,i_sw,i_oom) = mean(diff(pks2(:))./diff(locs2(:)));
                %-- oscillatopn time
                oscilT1_MAT(i_p,i_sw,i_oom) = mean(diff(locs1));
                oscilT2_MAT(i_p,i_sw,i_oom) = mean(diff(locs2));
                %-- amplitude
                Amp1_MAT(i_p,i_sw,i_oom) = peak2peak(phediDiff(timeVec>=(max(timeVec)-1.5*oscilT1),1));
                Amp2_MAT(i_p,i_sw,i_oom) = peak2peak(phediDiff(timeVec>=(max(timeVec)-1.5*oscilT1),2));
                
                pause(0.8);
            end
        end
    end
end

%% plot
if 0
    [S_W_grid, P_grod] = meshgrid(S_W,P);
    Plin = linspace(min(S_W),max(P),33);
    SWlin = linspace(min(S_W),max(P),33);
    [X,Y] = meshgrid(Plin,SWlin);
    
    figure;
    driftVel1_MAT_intrp = interp2(S_W_grid, P_grod,driftVel1_MAT(:,:,1),X,Y);
    driftVel2_MAT_intrp = interp2(S_W_grid, P_grod,driftVel1_MAT(:,:,2),X,Y);
    subplot(1,2,1); surf(S_W_grid,P_grod,driftVel1_MAT_intrp); xlabel('sterp width'); ylabel('interpolant P'); zlabel('drift velocity');
    subplot(1,2,2); surf(S_W_grid,P_grod,driftVel2_MAT_intrp); xlabel('sterp width'); ylabel('interpolant P'); zlabel('drift velocity');
    
    figure;
    oscilT1_MAT_intrp = interp2(S_W_grid, P_grod,oscilT1_MAT(:,:,1),X,Y);
    oscilT2_MAT_intrp = interp2(S_W_grid, P_grod,oscilT2_MAT(:,:,2),X,Y);
    subplot(1,2,1); surf(S_W_grid,P_grod,oscilT1_MAT_intrp(:,:,1)); xlabel('sterp width'); ylabel('interpolant P'); zlabel('scilatn T');
    subplot(1,2,2); surf(S_W_grid,P_grod,oscilT2_MAT(:,:,2)); xlabel('sterp width'); ylabel('interpolant P'); zlabel('scilatn T');
    figure;
    subplot(1,2,1); surf(S_W_grid,P_grod,Amp1_MAT(:,:,1)); xlabel('sterp width'); ylabel('interpolant P'); zlabel('Amp1_MAT');
    subplot(1,2,2); surf(S_W_grid,P_grod,Amp1_MAT(:,:,2)); xlabel('sterp width'); ylabel('interpolant P'); zlabel('Amp1_MAT');
end
%% save
folder = 'C:\Users\NeriB\Google Drive\JAY_lab\PROGRESS\phedi-simulation';
for iFig = 1:length(FigList)
    FigHandle = NewFigList(iFig);
    FigName   = ['P',num2str(P_names(iFig)),'S_W',num2str(S_W_names(iFig)),'OOM',num2str(numOOM_names(iFig))];
    FigName2 = strrep(FigName,'_',[]);
    FigName3 = strrep(FigName2,'.','_');
    fullName = fullfile(folder, [FigName3, '.fig']);
    savefig(FigHandle, fullName);
    export_fig(FigHandle,fullName(1:end-4),'-png','-r600')
end