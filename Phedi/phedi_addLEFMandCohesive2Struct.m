function NewDataStruct = phedi_addLEFMandCohesive2Struct(DataStruct,Param)

%% load data
if nargin<2
    Param = DataStruct.Param;
end
sgDataStruct = DataStruct.SgData;
BigPicRowOverTimeStruct = DataStruct.BigPicRotStruct;

%% turn off warning for this section
warning('off');

%% solution at strain gages
disp('Calculating Gamma...');
%--- find the gamma to fit to
relevantSGs = find(Param.scratchRegionMeters(1)<sgDataStruct.x_sg*1e-3 & sgDataStruct.x_sg*1e-3<Param.scratchRegionMeters(2));
relevantSG = round(mean(relevantSGs));  % take only on SG
% solAtSG = gamma_findGammaFit(sgDataStruct,BigPicRowOverTimeStruct, relevantSG);
solAtSG = gamma_findGammaFit_ByAmpUxx(sgDataStruct,BigPicRowOverTimeStruct, relevantSG);
solAtSG.SG_calc = relevantSG;
disp(['Gamma=',num2str(solAtSG.Gamma)]);

%% calc the cohesive zone model
disp('calculating cohesive zone model on interface...');
v = DataStruct.PhediData.Cf/Param.Cr;
% CohesiveModelStruct=CrackSolutionGeneralCohesive_giveParameters(v,Param.CohsvModl_h,solAtSG.Gamma);
% CohesiveModelStruct = CrackSolutionGeneralCohesive_InsertVariables(v,Param.CohsvModl_h,solAtSG.Gamma,'L0',3e-3,'Ohnaka89');
CohesiveModelStruct = [];
CohesiveModelStruct.h = Param.CohsvModl_h;
CohesiveModelStruct.v = v;
CohesiveModelStruct.Cr = Param.Cr;

%% calculate LEFM on interface
disp('calculating LEFM sol. on interface...');
interfaceEpsilon = 1e-9;
X =  -(-0.5:1e-6:0.5);
solAtInter = CrackSolutionForh_GammaChange(v, -0.5, interfaceEpsilon, solAtSG.Gamma, X);
% [power,sqrtPref] = Crack_UxFromInterfaceSolution(solAtInter);
% XX = -solAtInter.x;
% Ux = sqrtPref*(X.^power);
% Ux(XX>=0) = 0;
% solAtInter.Ux = Ux;
% solAtInter.x_4Ux = XX;

%% update dataStruct
NewDataStruct = DataStruct;
NewDataStruct.CohesiveModelStruct = CohesiveModelStruct;
NewDataStruct.solAtSG = solAtSG;
NewDataStruct.solAtInter = solAtInter;

%% turn back on warning messages
warning('on');


