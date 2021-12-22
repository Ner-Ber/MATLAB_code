function sol = gamma_findGammaFit_wShift(SgData,BigPicRotStruct, varargin)
% fitStruct = gamma_findGammaFit(SgData,BigPicRotStruct,sg_indexes,strainComp,Nsmth, xShifts)
%
%gamma_findGammaFit will find the Gamma (energy) of a measurement by
%least-squares fitting. done using only a single component of the strain
%tensor.

%% defaults and ata
% BigPicRotStruct = DataStruct.BigPicRotStruct;
% SgData = DataStruct.SgData;

[sg_indexes,strainComp,Nsmth, xShifts] = setDefaults4function(varargin,8,'Uxx',5, linspace(-0.02,0.02,300));

Cr = 1255;
%% fit
gamma_vec = 1:0.2:15;
x_lim = 0.03;

all_x = cell(length(gamma_vec),length(xShifts));
all_Strain = cell(length(gamma_vec),length(xShifts));
all_rms = nan(length(gamma_vec),length(xShifts));

%--- strain from sg:
if strcmp(strainComp,'Uxx')
    strain_vec = (SgData.Uxx(:,sg_indexes)-mean(SgData.Uxx(1:100,sg_indexes)))*1e-3;
elseif strcmp(strainComp,'Uxy')
    strain_vec = (SgData.Uxy(:,sg_indexes)-mean(SgData.Uxy(1:100,sg_indexes)))*1e-3;
elseif strcmp(strainComp,'Uyy')
    strain_vec = (SgData.Uyy(:,sg_indexes)-mean(SgData.Uyy(1:100,sg_indexes)))*1e-3;
end

%--- find x_min_x_tip:
[~,IdxMin] = min(abs(BigPicRotStruct.x-SgData.x_sg(sg_indexes)/1000));
SG_t_tip = BigPicRotStruct.frontTime_interp(IdxMin)/BigPicRotStruct.fps;
t_mins_t_tip_vec = (SgData.t./1000)-SG_t_tip;
x_mins_x_tip = -t_mins_t_tip_vec*SgData.v_sg(sg_indexes);

%--- if Uxy reference to other side
if strcmp(strainComp,'Uxy')
    reference_logical = (x_mins_x_tip<=-0.04) & (x_mins_x_tip>=(-0.06));
    strain_vec = (SgData.Uxy(:,sg_indexes)-mean(SgData.Uxy(reference_logical,sg_indexes)))*1e-3;
end

strain_smth = smooth(strain_vec,Nsmth);
%--- get current solution and iterate on gammas
v = SgData.v_sg(sg_indexes);
h = SgData.y_sg(sg_indexes)*1e-3;

for xs = 1:length(xShifts)
    thisX = x_mins_x_tip+xShifts(xs);
    relevant_logical = (thisX<=x_lim) & (thisX>=(-x_lim));
    relevantX = thisX(relevant_logical);
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
        
        all_Strain{g,xs} = solStrain;
        all_x{g,xs} = sol.x;
        
        A_S = abs(solStrain);
        A_S_norm = A_S-min(A_S);
        weights = (A_S_norm./max(A_S_norm))*0.8+0.2;    % weights for error calculation
        
        all_rms(g,xs) = rms(weights.*(strain_smth(relevant_logical)-solStrain));
%         all_rms(g,xs) = r_square(strain_smth(relevant_logical),solStrain);
    end
end

%--- find minimal rms:
[rms_mins, Idx] = min(all_rms(:));
% [rms_mins, Idx] = max(all_rms(:));
[g_idx,shift_idx] = ind2sub(size(all_rms),Idx);
%% build fit struct
sol = CrackSolutionForh_GammaChange(v/Cr,-0.5,h,gamma_vec(g_idx),-0.1:1e-4:0.1);
sol.rms_mins = rms_mins;
sol.shift_idx = shift_idx;
end

