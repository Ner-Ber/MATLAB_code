%% X_pk_LEFM simulation

[Cd, Cs, Cr]=CrackSolutionMaterialProperties;
Gamma = 3.5;
Cf_vec = linspace(1^2.8,1255^2.8,100).^(1/2.8);
L = 3.5e-3;
model = 'exp';
x_vec = [-2e-3:1e-4:(-8e-4-5e-6),-8e-4:5e-6:2e-4]; %-0.03:1e-5:0.01;
L_estimate_model = CohesiveSolutionGetL(Cf_vec./Cr,Gamma,2e6)/10;
maxCohesive = nan(size(Cf_vec));
X_pkLEFM = nan(size(Cf_vec));

for c = 1:length(Cf_vec)
    sol_cohesive = CrackSolutionGeneralCohesive_InsertVariables_Neri(Cf_vec(c)/Cr,1e-9,Gamma,'L',3.5e-3/10,model,[],x_vec);
    sol_LEFM = CrackSolutionForh_GammaChange(Cf_vec(c)/Cr , -0.5 , 1e-9 , Gamma, x_vec);
    
    maxCohesive(c) = max(sol_cohesive.vx(:));
    [~,I] = min(abs(sol_LEFM.vx-maxCohesive(c)));
    X_pkLEFM(c) = sol_LEFM.x(I);
    
end