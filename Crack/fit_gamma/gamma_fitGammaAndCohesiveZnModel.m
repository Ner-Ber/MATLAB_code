function [solAtSG,CohesiveModelStruct] = gamma_fitGammaAndCohesiveZnModel(DataStruct,varargin)

%% defaults
doPlot = setDefaults4function(varargin,0);

Cr = 1255;
CohsvModl_h = 1e-7;
if ~isfield(DataStruct.ExperimentData,'scratchRegionMeters')
    scratchRegionMeters = [0.07 0.09];
else
    scratchRegionMeters = DataStruct.ExperimentData.scratchRegionMeters;
end

%%  find the gamma to fit to
relevantSGs = find(scratchRegionMeters(1)<DataStruct.SgData.x_sg*1e-3 & DataStruct.SgData.x_sg*1e-3<scratchRegionMeters(2));
relevantSG = round(mean(relevantSGs));  % take only on SG
solAtSG = gamma_findGammaFit(DataStruct.SgData,DataStruct.BigPicRotStruct, relevantSG);
solAtSG.SG_calc = relevantSG;
disp(['Gamma=',num2str(solAtSG.Gamma)]);


%% calc the cohesive zone model
disp('calculating cohesive zone model on interface...');
% [~, ~, Cr, ~ , ~, ~, ~,~,~,~]=CrackSolutionMaterialProperties;

v = DataStruct.BigPicRotStruct.AvgCf_scratches/Cr;
CohesiveModelStruct=CrackSolutionGeneralCohesive_giveParameters(v,CohsvModl_h,solAtSG.Gamma);
% CohesiveModelStruct=CrackSolutionGeneralCohesive(v,h);
CohesiveModelStruct.h = CohsvModl_h;
CohesiveModelStruct.v = v;
CohesiveModelStruct.Cr = Cr;


%% plot
if doPlot
    fakeStruct = struct;
    fakeStruct.SgData = DataStruct.SgData;
    fakeStruct.solAtSG = solAtSG;
    gamma_plotGammaFit(fakeStruct);
end

end