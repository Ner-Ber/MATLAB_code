function sol = gamma_findGammaFit_ByAmpUxx(SgData,BigPicRotStruct, varargin)
    % fitStruct = gamma_findGammaFit(SgData,BigPicRotStruct,sg_indexes,x_lim)
    %
    %gamma_findGammaFit will find the Gamma (energy) of a measurement by
    %least-squares fitting. done using only a single component of the strain
    %tensor.
    
    %% defaults and ata
    % BigPicRotStruct = DataStruct.BigPicRotStruct;
    % SgData = DataStruct.SgData;
    
    [sg_indexes,x_lim,strainComp,Nsmth] = setDefaults4function(varargin,8,0.04,'Uxx',5);
    
    [Cd, Cs, Cr] = CrackSolutionMaterialProperties;
    %% fit
    gamma_vec = 1:0.2:15;
    %--- strain from sg:
    % strain_vec = (SgData.Uxx(:,sg_indexes)-mean(SgData.Uxx(1:100,sg_indexes)))*1e-3;
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
    x_mins_x_tip = signal_ChangeTime2Space_mapping(t_mins_t_tip_vec*1e-3, SgData.x_sg(sg_indexes)/1000, BigPicRotStruct);
    % if isfield(SgData,'x_mins_x_tips')
    %     x_mins_x_tip = SgData.x_mins_x_tips(:,sg_indexes);
    % else
    %     [~,IdxMin] = min(abs(BigPicRotStruct.x-SgData.x_sg(sg_indexes)/1000));
    %     SG_t_tip = BigPicRotStruct.frontTime_interp(IdxMin)/BigPicRotStruct.fps;
    %     t_mins_t_tip_vec = (SgData.t./1000)-SG_t_tip;
    %     x_mins_x_tip = signal_ChangeTime2Space_mapping(t_mins_t_tip_vec*1e-3, SgData.x_sg(sg_indexes)/1000, BigPicRotStruct);
    % %     x_mins_x_tip = -t_mins_t_tip_vec*SgData.v_sg(sg_indexes);
    % end
    
    
    %--- select region where to check for high amp
    relevant_logical = (x_mins_x_tip<=x_lim) & (x_mins_x_tip>=(-x_lim));
    X_reduced = x_mins_x_tip(relevant_logical);
    
    %--- get location of maximal amplitude
    % strain_smth = supsmu(X_reduced,strain_vec(relevant_logical));
    strain_smth = smooth(strain_vec(relevant_logical),Nsmth);
    [maxG] = max(abs(strain_smth));
    
    %--- get current solution and iterate on gammas
    v = SgData.v_sg(sg_indexes);
    if v==inf || v==-inf
        v = mean(SgData.v_sg([sg_indexes-1 sg_indexes+1]));
    end
    h = SgData.y_sg(sg_indexes)*1e-3;
    maxG_vec = nan(1,length(gamma_vec));
    for g = 1:length(gamma_vec)
        sol = CrackSolutionForh_GammaChange(v/Cr,-0.5,h,gamma_vec(g),X_reduced);
        maxG_vec(g) = max(abs(sol.Uxx));
    end
    
    if maxG<=max(maxG_vec) && maxG>=min(maxG_vec)
        [FoundGamma,~] = intersections(gamma_vec,maxG_vec,[-100 100],[1 1]*maxG);
    else
        FoundGamma = interp1(maxG_vec,gamma_vec,maxG,'linear','extrap');
    end
    
    %% build fit struct
    sol = CrackSolutionForh_GammaChange(v/Cr,-0.5,h,FoundGamma,-0.1:1e-4:0.1);
end

