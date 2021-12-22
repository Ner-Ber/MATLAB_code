%% correlate between shear and particle velocity on interface by cohesive zone model


[Cd, Cs, Cr, nu, ro, E, mu, ~, PlaneStrain, ~, ~]=CrackSolutionMaterialProperties;
Model = 'exp';
% L_vec = (0.4:0.5:5)*1e-3;
% Gamma_vec = 0.8:0.5:5;
% v_vec = 50:150:1250;

L_vec = (5)*1e-3;
Gamma_vec = 2.5;
v_vec = [1:10:1000, 1002:2:1255];

velPeakMat = nan(length(L_vec),length(v_vec),length(Gamma_vec));
velLockMat = nan(length(L_vec),length(v_vec),length(Gamma_vec));
tauPeakMat = nan(length(L_vec),length(v_vec),length(Gamma_vec));

for l = 1:length(L_vec)
    for v = 1:length(v_vec)
        for g = 1:length(Gamma_vec)
            [X_peak, Y_peak] = cohesive_findPeakOfVelByIterations(v_vec(v)/Cr,1e-8,Gamma_vec(g),'L',L_vec(l),Model);
            tau_p = CohesiveSolutionGetTau_p(v_vec(v)/Cr,Gamma_vec(g),L_vec(l),Model);
            
            velPeakMat(l,v,g) = Y_peak;
            velLockMat(l,v,g) = X_peak;
            tauPeakMat(l,v,g) = tau_p;
        end
    end
end

[L,V,G] = meshgrid(L_vec,v_vec,Gamma_vec);
ratioVals = velPeakMat./tauPeakMat;
% figure; scatter3(L(:),V(:),G(:),[],ratioVals(:),'filled');
% xlabel('L [m]');
% ylabel('C_f [m/s]');
% zlabel('\Gamma [J/m^2]');


figure; hold on;
for l = 1:length(L_vec)
plot(Gamma_vec,reshape(ratioVals(l,5,:),[],1),'DisplayName',['L=',num2str(l),'  C_f=',num2str(v_vec(5))]);
end

figure; hold on;
for g = 1:length(Gamma_vec)
plot(v_vec,reshape(ratioVals(5,:,g),[],1),'DisplayName',['L=',num2str(L_vec(5)),'  \Gamma=',num2str(Gamma_vec(g))]);
end

figure; hold on;
for v = 1:length(v_vec)
plot(L_vec,reshape(ratioVals(:,v,5),[],1),'DisplayName',['C_f=',num2str(v_vec(v)),'  \Gamma=',num2str(Gamma_vec(5))]);
end
