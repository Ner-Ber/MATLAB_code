%% analyze errors of phedi algorithm from simulation
%% to create the data i used:
% PhediLocationPixSprMat = [];
% shiftVectorMat = [];
% for Vslip=VslipVec
%     follow_phedi_simulation
%     pause(0.1);
%     PhediLocationPixSprMat = cat(3,PhediLocationPixSprMat,PhediLocationPix);
%     shiftVectorMat = cat(2,shiftVectorMat,shiftVector(firstFrame4Phedi:lastFrame4Phedi)'*dx);
% end

%% load data from simulation with details: method=minimal RMS     N-OOM=3      interfrnce dist.=0.002
load('C:\Users\NeriB\Google Drive\JAY_lab\FIGURES\algorithm and errors\data4errs\VelShiftAndMeasure.mat');
movementStarts = 32;
%% prepare relative difference measre-real
relativePhediDisplacement = bsxfun(@minus,PhediLocationPixSprMat,PhediLocationPixSprMat(1,:,:));
comapredPhediMovement = bsxfun(@minus,relativePhediDisplacement,permute(shiftVectorMat,[1 3 2]));
%% fit linear to each traj
L = movementStarts:size(PhediLocationPixSprMat,1);
P = zeros(2,size(PhediLocationPixSprMat,2),size(PhediLocationPixSprMat,3));
oscilationRemainMat = zeros(length(L),size(PhediLocationPixSprMat,2),size(PhediLocationPixSprMat,3));

for phediIdx = 1:size(PhediLocationPixSprMat,2)
    for Vidx = 1:length(VslipVec)
        p = polyfit(shiftVectorMat(L,Vidx),comapredPhediMovement(L,phediIdx,Vidx),1);
        P(:,phediIdx,Vidx) = p(:);
        oscilationRemain = comapredPhediMovement(L,phediIdx,Vidx) - polyval(p,shiftVectorMat(L,Vidx));
        oscilationRemainMat(:,phediIdx,Vidx) = oscilationRemain(:);
    end
end

%% plot drift
figure; hold on;
simStr = {'*-','.-'};
Vcolor = MyVaryColor(length(VslipVec));
for phediIdx = 1:size(PhediLocationPixSprMat,2)
    for Vidx = 1:length(VslipVec)
        plot(shiftVectorMat(L,Vidx),comapredPhediMovement(L,phediIdx,Vidx),...
            simStr{phediIdx},...
            'Color',Vcolor(Vidx,:),...
            'DisplayName',['phedi=',num2str(phediIdx),'  V_{slip}=',num2str(VslipVec(Vidx))]);
    end
end
ylabel('\Delta=(measured shift)-(real shift) [pix]');
xlabel('real shift [pix]')
title('accumulated drift')


%% plot oscillation remainder
figure; hold on;
simStr = {'*-','.-'};
Vcolor = MyVaryColor(length(VslipVec));
for phediIdx = 1:size(PhediLocationPixSprMat,2)
    for Vidx = 1:length(VslipVec)
        plot(shiftVectorMat(L,Vidx),oscilationRemainMat(:,phediIdx,Vidx),...
            simStr{phediIdx},...
            'Color',Vcolor(Vidx,:),...
            'DisplayName',['phedi=',num2str(phediIdx),'  V_{slip}=',num2str(VslipVec(Vidx))]);
    end
end
ylabel('local \Delta [pix]');
xlabel('real shift [pix]')
title('oscillation remainer')

%% compute and display errors
%-- oscillation errors
oscil = peak2peak(oscilationRemainMat(1:100,:,:),1);
oscil = squeeze(oscil);
figure; plot(VslipVec,oscil');
legend({'phedi1','phedi2'});
ylabel('oscilaltion p2p [pix]');
xlabel('Vslip [m/s]');

%-- drift errors
Delta2shift = P(1,:,:);
Delta2shift = squeeze(Delta2shift);
figure; plot(VslipVec,Delta2shift');
legend({'phedi1','phedi2'});
ylabel('\Delta/shift [pix/pix]');
xlabel('Vslip [m/s]');

