function cohesive_fitLengthByModel(DataStruct)

%% prepare data for Ohanaka89 model
l_vec = (0.1:0.4:5)*1e-3;
dx = linspace(5e-6,1e-4,length(l_vec));
x_range = [-0.02 0.01];
% x = max(x_range):-abs(dx):min(x_range);
all_Ux = cell(1,length(l_vec));
all_Uxx = cell(1,length(l_vec));
for l=1:length(l_vec)
    thisSol = CrackSolutionGeneralCohesive_InsertVariables(785/1255,2e-8,PhediStructCell{6}.solAtSG.Gamma,'L0',l_vec(l),'Ohnaka89',[],x_range,dx(l));
    all_Ux{l} = thisSol.Ux(:);
    all_Uxx{l} = thisSol.Uxx(:);
end

%%
relevantPhedis = find(DataStruct.PhediData.slopeIncline==1);
[T_mean,PhediLoc_mean,PhediLoc_max,PhediLoc_min] = phedi_averagePhedi(DataStruct, 'vel', relevantPhedis);
X_mean = signal_ChangeTime2Space_mapping(T_mean, DataStruct.SgData.x_sg(8)*1e-3, DataStruct.BigPicRotStruct);


X_mean,PhediLoc_mean

