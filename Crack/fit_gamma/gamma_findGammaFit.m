function sol = gamma_findGammaFit(SgData,BigPicRotStruct, varargin)
% fitStruct = gamma_findGammaFit(DataStruct, sg_indexes, Nsmth)
%
%gamma_findGammaFit will find the Gamma (energy) of a measurement by
%least-squares fitting. done using only a single component of the strain
%tensor.

%% defaults and ata
% BigPicRotStruct = DataStruct.BigPicRotStruct;
% SgData = DataStruct.SgData;

[sg_indexes,strainComp,Nsmth] = setDefaults4function(varargin,8,'Uxx',15);

Cr = 1255;
%% fit
gamma_vec = 1:0.2:15;
x_lim = 0.05;

all_x = cell(length(gamma_vec),1);
all_Strain = cell(length(gamma_vec),1);
all_rms = nan(length(gamma_vec),1);

%--- strain from sg:
if strcmp(strainComp,'Uxx')
    strain_vec = (SgData.Uxx(:,sg_indexes)-mean(SgData.Uxx(1:100,sg_indexes)))*1e-3;
elseif strcmp(strainComp,'Uxy')
    strain_vec = (SgData.Uxy(:,sg_indexes)-mean(SgData.Uxy(1:100,sg_indexes)))*1e-3;
elseif strcmp(strainComp,'Uyy')
    strain_vec = (SgData.Uyy(:,sg_indexes)-mean(SgData.Uyy(1:100,sg_indexes)))*1e-3;
end

%--- find x_min_x_tip:
if isfield(SgData,'x_mins_x_tips')
    x_mins_x_tip = SgData.x_mins_x_tips(:,sg_indexes);
else
    [~,IdxMin] = min(abs(BigPicRotStruct.x-SgData.x_sg(sg_indexes)/1000));
    SG_t_tip = BigPicRotStruct.frontTime_interp(IdxMin)/BigPicRotStruct.fps;
    t_mins_t_tip_vec = (SgData.t./1000)-SG_t_tip;
    x_mins_x_tip = -t_mins_t_tip_vec*SgData.v_sg(sg_indexes);
end

%--- if Uxy reference to other side
if strcmp(strainComp,'Uxy')
    reference_logical = (x_mins_x_tip<=-0.04) & (x_mins_x_tip>=(-0.06));
    strain_vec = (SgData.Uxy(:,sg_indexes)-mean(SgData.Uxy(reference_logical,sg_indexes)))*1e-3;
end

strain_smth = smooth(strain_vec,Nsmth);
%--- get current solution and iterate on gammas
v = SgData.v_sg(sg_indexes);
h = SgData.y_sg(sg_indexes)*1e-3;
relevant_logical = (x_mins_x_tip<=x_lim) & (x_mins_x_tip>=(-x_lim));
relevantX = x_mins_x_tip(relevant_logical);

for g = 1:length(gamma_vec)
    sol = CrackSolutionForh_GammaChange(v/Cr,-0.5,h,gamma_vec(g),relevantX);
    if strcmp(strainComp,'Uxx')
        solStrain = sol.Uxx;
    elseif strcmp(strainComp,'Uxy')
        solStrainTMP = sol.Uxy;
        referenceSOL_logical = (sol.x<=-0.04) & (sol.x>=(-0.06));
        solStrain = solStrainTMP(:)-mean(solStrainTMP(referenceSOL_logical));
    elseif strcmp(strainComp,'Uyy')
        solStrain = sol.Uyy;
    end
    
    all_Strain{g} = solStrain;
    all_x{g} = sol.x;
    
    all_rms(g) = rms(strain_smth(relevant_logical)-solStrain);
end

%--- find minimal rms:
[rms_mins, g_idx] = min(all_rms);
%% build fit struct
sol = CrackSolutionForh_GammaChange(v/Cr,-0.5,h,gamma_vec(g_idx),-0.1:1e-4:0.1);
sol.rms_mins = rms_mins;
end

