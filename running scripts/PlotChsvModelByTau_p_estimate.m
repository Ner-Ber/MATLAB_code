%% fix phediStructCell
relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
relevantEvents = relevantEvents(relevantEvents~=1)
for j=1:length(relevantEvents)
    PhediStructCell{relevantEvents(j)}.PhediDataSG = PhediStructCell{relevantEvents(j)}.PhediData;
    PhediStructCell{relevantEvents(j)}.PhediData = PhediStructCell{relevantEvents(j)}.PhediDataSimpSmth;
end
%%
[Cd, Cs, Cr, nu, ro, E, mu, Gamma, PlaneStrain, tau_p, Xc0]=CrackSolutionMaterialProperties;

%% get relevant L and plot
A_avgDist = 0.01;
A_ratio = 0.02;
model = 'exp';

L_vec = nan(1,length(relevantEvents));
tau_p = nan(1,length(relevantEvents));
for i=1:length(relevantEvents)
    try
        yeildSall = yield_findYieldStressStrain(PhediStructCell{relevantEvents(i)},A_avgDist);
        yieldSxy2nd = yeildSall.yieldSxy2nd;
        Cf = PhediStructCell{relevantEvents(i)}.PhediData.Cf;
        Gamma = PhediStructCell{relevantEvents(i)}.solAtSG.Gamma;
        tau_p(i) = yieldSxy2nd./A_ratio;
        L_vec(i) = 10*CohesiveSolutionGetL(Cf./Cr,Gamma,tau_p(i),model);
        
        sol = CrackSolutionGeneralCohesive_InsertVariables_Neri(Cf/Cr,1e-9,Gamma,'L',L_vec(i),model);
        Movie_phedi_plotWithSg(PhediStructCell{relevantEvents(i)},'00',{'avgPhediVel'});
        hold on;
        plot(sol.x, sol.vx);
    catch
    end
end


