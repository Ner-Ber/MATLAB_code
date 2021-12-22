%% contact area, v_peak and residual stress

D = dir_names;
% k = strfind(D,'.mat');
k1 = strfind(D,'.mat')
k2 = strfind(D,'allPhedi')
% F = find(~cellfun(@isempty,k))
F = find(~cellfun(@isempty,k1) & ~cellfun(@isempty,k2))
%% parameters
residual_locs = -0.1:0.01:-0.01;
[~,~, Cr]=CrackSolutionMaterialProperties;
myStruct = struct([]);
A_avgDist = 0.01;
for j=1:length(F)
    %% iterate on events
    try
        disp(['loading ',D{F(j)}]);
        load(D{F(j)});
        
        relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
        relevantEvents = relevantEvents(relevantEvents~=1);
        
        Cf_vec = zeros(length(relevantEvents),1);
        Gamma_vec = zeros(length(relevantEvents),1);
        DA_vec_tmp = zeros(length(residual_locs),length(relevantEvents));
        DA_vec_spt = zeros(length(residual_locs),length(relevantEvents));
        DA_vec_sml = zeros(length(residual_locs),length(relevantEvents));
        Ar_vec_tmp = zeros(length(residual_locs),length(relevantEvents));
        Ar_vec_spt = zeros(length(residual_locs),length(relevantEvents));
        Ar_vec_sml = zeros(length(residual_locs),length(relevantEvents));
        peakVel_vecAll = zeros(length(relevantEvents),1);
        peakVel_vecPos = zeros(length(relevantEvents),1);
        peakVel_vecNeg = zeros(length(relevantEvents),1);
        
        VarVPeak_vecAll = zeros(length(relevantEvents),1);
        VarVPeak_vecPos = zeros(length(relevantEvents),1);
        VarVPeak_vecNeg = zeros(length(relevantEvents),1);
        
        peakLoc_vecAll = zeros(length(relevantEvents),1);
        peakLoc_vecPos = zeros(length(relevantEvents),1);
        peakLoc_vecNeg = zeros(length(relevantEvents),1);
        
        
        FWHM_vecAll = zeros(length(relevantEvents),1);
        FWHM_vecPos = zeros(length(relevantEvents),1);
        FWHM_vecNeg = zeros(length(relevantEvents),1);
        
        yieldSxy1stAll = zeros(length(relevantEvents),1);
        yieldSxy1stPos = zeros(length(relevantEvents),1);
        yieldSxy1stNeg = zeros(length(relevantEvents),1);
        yieldSxy2ndAll = zeros(length(relevantEvents),1);
        yieldSxy2ndPos = zeros(length(relevantEvents),1);
        yieldSxy2ndNeg = zeros(length(relevantEvents),1);
        
        rise1XAll = zeros(length(relevantEvents),1);
        rise1XPos = zeros(length(relevantEvents),1);
        rise1XNeg = zeros(length(relevantEvents),1);
        
        rise2XAll = zeros(length(relevantEvents),1);
        rise2XPos = zeros(length(relevantEvents),1);
        rise2XNeg = zeros(length(relevantEvents),1);
        rise2VAll = zeros(length(relevantEvents),1);
        rise2VPos = zeros(length(relevantEvents),1);
        rise2VNeg = zeros(length(relevantEvents),1);
        
        Vresiduals_matAll = zeros(length(residual_locs),length(relevantEvents));
        Vresiduals_matPos = zeros(length(residual_locs),length(relevantEvents));
        Vresiduals_matNeg = zeros(length(residual_locs),length(relevantEvents));
        
        shifts_cell = cell(length(relevantEvents),1);
        RMS_cell = cell(length(relevantEvents),1);
        
        tauResiduals_mat = zeros(length(residual_locs),length(relevantEvents));
        
        [maxXi_loc_vec,maxYi_loc_vec,maxXi_vel_vec,maxYi_vel_vec] = ...
            deal(zeros(length(relevantEvents),1));
        
        for iEv = 1:length(relevantEvents)
            try
                
                PhediStructCell{relevantEvents(iEv)}.PhediDataSG = PhediStructCell{relevantEvents(iEv)}.PhediData;
                PhediStructCell{relevantEvents(iEv)}.PhediData = PhediStructCell{relevantEvents(iEv)}.PhediDataSimpSmth;
                
                BigPicRotStruct = PhediStructCell{relevantEvents(iEv)}.BigPicRotStruct;
                PhediData = PhediStructCell{relevantEvents(iEv)}.PhediData;
                slopeIncline = mean(PhediStructCell{relevantEvents(iEv)}.PhediData.slopeIncline,1);
                %% get DA
                %--- get contact area
                x_minXtip_4spatial = BigPicRotStruct.x-PhediData.PhotoLocation;
                frontTime = BigPicRotStruct.frontTime_interp./BigPicRotStruct.fps;
                [~,I1] = min(abs(x_minXtip_4spatial));
                t_minTtip_4tmprl = BigPicRotStruct.t-frontTime(I1);
                [~,I2] = min(abs(t_minTtip_4tmprl));
                A_spatial = BigPicRotStruct.DataMatNorm(I2,:);
                A_temporal_old = BigPicRotStruct.DataMatNorm(:,I1);
                x_minXtip_4tmprl_old = signal_ChangeTime2Space_mapping(t_minTtip_4tmprl,PhediData.PhotoLocation,BigPicRotStruct);
                cropLogic = x_minXtip_4tmprl_old>=-0.1 & x_minXtip_4tmprl_old<=0.1;
                x_minXtip_4tmprl_old = x_minXtip_4tmprl_old(cropLogic);
                A_temporal_old = A_temporal_old(cropLogic);
                x_minXtip_4tmprl = max(x_minXtip_4tmprl_old):-1e-4:min(x_minXtip_4tmprl_old);
                A_temporal = interp1(x_minXtip_4tmprl_old,A_temporal_old,x_minXtip_4tmprl);
                
                ROTsum = sum(PhediStructCell{relevantEvents(iEv)}.AsperityData.RowOverTimePreBlur,2);
                A_smallPic_old = ROTsum./mean(ROTsum(1:1000));
                t_minTtips_small = PhediStructCell{relevantEvents(iEv)}.AsperityData.timeCount - mean(PhediData.t_tips);
                x_minXtip_4small_old = signal_ChangeTime2Space_mapping(t_minTtips_small,PhediData.PhotoLocation,BigPicRotStruct);
                cropLogic = x_minXtip_4small_old>=-0.1 & x_minXtip_4small_old<=0.1;
                x_minXtip_4small_old = x_minXtip_4small_old(cropLogic);
                A_smallPic_old = A_smallPic_old(cropLogic);
                x_minXtip_4small = max(x_minXtip_4small_old):-1e-4:min(x_minXtip_4small_old);
                A_smallPic = interp1(x_minXtip_4small_old,A_smallPic_old,x_minXtip_4small);
                
                %--- smooth contact area
                dxTmprl =abs(diff(x_minXtip_4tmprl([1,2])));
                smthPixTmprl = A_avgDist./dxTmprl;
                A_temporalSmth = smooth(A_temporal, smthPixTmprl);
                dxSptl = abs(diff(x_minXtip_4spatial([1,2])));
                smthPixSptal = A_avgDist./dxSptl;
                A_spatialSmth= smooth(A_spatial, smthPixSptal);
                dxSmall = abs(diff(x_minXtip_4small([1,2])));
                smthPixSmall = A_avgDist./dxSmall;
                A_SmallSmth= smooth(A_smallPic, smthPixSmall);
                
                
                %--- get DA
                [~,Itmp] = min(abs(bsxfun(@minus,x_minXtip_4tmprl(:),residual_locs(:)')));
                ArTmprl = A_temporalSmth(Itmp);
                [~,Ispt] = min(abs(bsxfun(@minus,x_minXtip_4spatial(:),residual_locs(:)')));
                ArSptl = A_spatialSmth(Ispt);
                [~,Isml] = min(abs(bsxfun(@minus,x_minXtip_4small(:),residual_locs(:)')));
                ArSmall = A_SmallSmth(Isml);
                
                DAtmprl = 1 - ArTmprl;
                DAsptl = 1 - ArSptl;
                DAsmall = 1 - ArSmall;
                
                %% get peak velocity from all phedis
                Sall = phedi_FWHMofMean(PhediStructCell{relevantEvents(iEv)},A_avgDist);
                VpeakXall = Sall.peakX;
                AvgPhediStructAll = Sall.AvgPhediStruct;
                %--- peak vel
                [~,I] = min(abs(AvgPhediStructAll.X_mean- VpeakXall));
                peakVelall = AvgPhediStructAll.Vel_mean_x(I);
                %-- var at peak
                PeakVelVarall = Sall.AvgPhediStruct.Vel_var_x(I);
                %--- velocity residulas
                [~,I] = min(abs(bsxfun(@minus,AvgPhediStructAll.X_mean(:),residual_locs(:)')));
                VresidualsAll = AvgPhediStructAll.Vel_mean_x(I);
                
                %-- yeilding point and kink
                yeildSall = yield_findYieldStressStrain(PhediStructCell{relevantEvents(iEv)},A_avgDist);
                
                %% get velocity properties from positive slope phedis
                posSlp = slopeIncline==+1;
                Spos = phedi_FWHMofMean(PhediStructCell{relevantEvents(iEv)},A_avgDist,posSlp);
                VpeakXpos = Spos.peakX;
                AvgPhediStructPos = Spos.AvgPhediStruct;
                %--- peak vel
                [~,I] = min(abs(AvgPhediStructPos.X_mean- VpeakXpos));
                peakVelPos = AvgPhediStructPos.Vel_mean_x(I);
                %-- var at peak
                PeakVelVarpos = Spos.AvgPhediStruct.Vel_var_x(I);
                %--- velocity residulas
                [~,I] = min(abs(bsxfun(@minus,AvgPhediStructPos.X_mean(:),residual_locs(:)')));
                VresidualsPos = AvgPhediStructPos.Vel_mean_x(I);
                
                %-- yeilding point and kink
                yeildSpos = yield_findYieldStressStrain(PhediStructCell{relevantEvents(iEv)},A_avgDist,posSlp);
                
                %% get velocity properties from negative slope phedis
                negSlp = slopeIncline==-1;
                Sneg = phedi_FWHMofMean(PhediStructCell{relevantEvents(iEv)},A_avgDist,negSlp);
                VpeakXneg = Sneg.peakX;
                AvgPhediStructNeg = Sneg.AvgPhediStruct;
                %--- peak vel
                [~,I] = min(abs(AvgPhediStructNeg.X_mean- VpeakXneg));
                peakVelNeg = AvgPhediStructNeg.Vel_mean_x(I);
                %-- var at peak
                PeakVelVarneg = Sneg.AvgPhediStruct.Vel_var_x(I);
                %--- velocity residulas
                [~,I] = min(abs(bsxfun(@minus,AvgPhediStructNeg.X_mean(:),residual_locs(:)')));
                VresidualsNeg = AvgPhediStructNeg.Vel_mean_x(I);
                
                %-- yeilding point and kink
                yeildSneg = yield_findYieldStressStrain(PhediStructCell{relevantEvents(iEv)},A_avgDist,negSlp);
                
                %% get tau residual
                X_MIN_xTIP_sgOLD = PhediStructCell{relevantEvents(iEv)}.SgData.x_mins_x_tips(:,8);
                cropLogic = X_MIN_xTIP_sgOLD>=-0.1 & X_MIN_xTIP_sgOLD<=0.1;
                X_MIN_xTIP_sgOLD = X_MIN_xTIP_sgOLD(cropLogic);
                Sxy_old = PhediStructCell{relevantEvents(iEv)}.SgData.Sxy(:,8);
                Sxy_old = Sxy_old(cropLogic);
                X_MIN_xTIP_sg = max(X_MIN_xTIP_sgOLD):-1e-4:min(X_MIN_xTIP_sgOLD);
                Sxy = interp1(X_MIN_xTIP_sgOLD,Sxy_old,X_MIN_xTIP_sg);
                smthL = A_avgDist/abs(diff(X_MIN_xTIP_sg([1,2])));
                SxySmth = smooth(Sxy,smthL);
                
                [~,I] = min(abs(bsxfun(@minus,X_MIN_xTIP_sg(:),residual_locs(:)')));
                tauResiduals = SxySmth(I);
                
                
                %% get delta_i (vertical shofts of Ux)
                [residualsCell, RMS_vec, shifts_vec] = phedi_UxBestFitVerticalShift(PhediStructCell{relevantEvents(iEv)});
                
                %% negative intersections with LEFM
                [maxXi_loc,maxYi_loc,maxXi_vel,maxYi_vel] =...
                    phedi_intersections_w_LEFM(PhediStructCell{relevantEvents(iEv)});
                
                

                %% keep stuff
                Cf_vec(iEv) = abs(PhediStructCell{relevantEvents(iEv)}.PhediData.Cf);
                Gamma_vec(iEv) = PhediStructCell{relevantEvents(iEv)}.solAtSG.Gamma;
                DA_vec_tmp(:,iEv) = DAtmprl;
                DA_vec_spt(:,iEv) = DAsptl;
                DA_vec_sml(:,iEv) = DAsmall;
                Ar_vec_tmp(:,iEv) = ArTmprl;
                Ar_vec_spt(:,iEv) = ArSptl;
                Ar_vec_sml(:,iEv) = ArSmall;
                peakVel_vecAll(iEv) = peakVelall;
                peakVel_vecPos(iEv) = peakVelPos;
                peakVel_vecNeg(iEv) = peakVelNeg;
                
                VarVPeak_vecAll(iEv) = PeakVelVarall;
                VarVPeak_vecPos(iEv) = PeakVelVarpos;
                VarVPeak_vecNeg(iEv) = PeakVelVarneg;
                
                peakLoc_vecAll(iEv) = VpeakXall;
                peakLoc_vecPos(iEv) = VpeakXpos;
                peakLoc_vecNeg(iEv) = VpeakXneg;
                FWHM_vecAll(iEv) = Sall.FWHM;
                FWHM_vecPos(iEv) = Spos.FWHM;
                FWHM_vecNeg(iEv) = Sneg.FWHM;
                yieldSxy1stAll(iEv) = yeildSall.yieldSxy1st;
                yieldSxy1stPos(iEv) = yeildSpos.yieldSxy1st;
                yieldSxy1stNeg(iEv) = yeildSneg.yieldSxy1st;
                yieldSxy2ndAll(iEv) = yeildSall.yieldSxy2nd;
                yieldSxy2ndPos(iEv) = yeildSpos.yieldSxy2nd;
                yieldSxy2ndNeg(iEv) = yeildSneg.yieldSxy2nd;
                
                rise1XAll(iEv) = yeildSall.risingVelPoint;
                rise1XPos(iEv) = yeildSpos.risingVelPoint;
                rise1XNeg(iEv) = yeildSneg.risingVelPoint;
                
                rise2XAll(iEv) = yeildSall.VelKinkX;
                rise2XPos(iEv) = yeildSpos.VelKinkX;
                rise2XNeg(iEv) = yeildSneg.VelKinkX;
                rise2VAll(iEv) = yeildSall.VelKinkY;
                rise2VPos(iEv) = yeildSpos.VelKinkY;
                rise2VNeg(iEv) = yeildSneg.VelKinkY;
                                
                Vresiduals_matAll(:,iEv) = VresidualsAll(:);
                Vresiduals_matPos(:,iEv) = VresidualsPos(:);
                Vresiduals_matNeg(:,iEv) = VresidualsNeg(:);
                
                tauResiduals_mat(:,iEv) = tauResiduals(:);
                
                shifts_cell{iEv} = shifts_vec;
                RMS_cell{iEv} = RMS_vec;
                
                [maxXi_loc_vec(iEv),maxYi_loc_vec(iEv),maxXi_vel_vec(iEv),maxYi_vel_vec(iEv)] = ...
                    deal(maxXi_loc,maxYi_loc,maxXi_vel,maxYi_vel);
            catch
                disp(['problem with some variables, Event=',num2str(relevantEvents(iEv))]);
            end
        end
        Name = [PhediStructCell{relevantEvents(iEv)}.ExperimentData.ExpDate,' ',PhediStructCell{relevantEvents(iEv)}.ExperimentData.ExpHour];
        
        myStruct(j).Cf_vec = Cf_vec;
        myStruct(j).Gamma = Gamma_vec;
        myStruct(j).DA_vec_tmp = DA_vec_tmp;
        myStruct(j).DA_vec_spt = DA_vec_spt;
        myStruct(j).DA_vec_sml = DA_vec_sml;
        myStruct(j).Ar_vec_tmp = Ar_vec_tmp;
        myStruct(j).Ar_vec_spt = Ar_vec_spt;
        myStruct(j).Ar_vec_sml = Ar_vec_sml;
        myStruct(j).peakVel_vecAll = peakVel_vecAll;
        myStruct(j).peakVel_vecPos = peakVel_vecPos;
        myStruct(j).peakVel_vecNeg = peakVel_vecNeg;
        myStruct(j).peakLoc_vecAll = peakLoc_vecAll;
        myStruct(j).peakLoc_vecPos = peakLoc_vecPos;
        myStruct(j).peakLoc_vecNeg = peakLoc_vecNeg;
        myStruct(j).VarVpeak_vecAll = VarVPeak_vecAll;
        myStruct(j).VarVpeak_vecPos = VarVPeak_vecPos;
        myStruct(j).VarVpeak_vecNeg = VarVPeak_vecNeg;
        myStruct(j).yieldSxy1stAll = yieldSxy1stAll;
        myStruct(j).yieldSxy1stPos = yieldSxy1stPos;
        myStruct(j).yieldSxy1stNeg = yieldSxy1stNeg;
        myStruct(j).yieldSxy2ndAll = yieldSxy2ndAll;
        myStruct(j).yieldSxy2ndPos = yieldSxy2ndPos;
        myStruct(j).yieldSxy2ndNeg = yieldSxy2ndNeg;
        myStruct(j).rise1XAll = rise1XAll;
        myStruct(j).rise1XPos = rise1XPos;
        myStruct(j).rise1XNeg = rise1XNeg;
        myStruct(j).rise2XAll = rise2XAll;
        myStruct(j).rise2XPos = rise2XPos;
        myStruct(j).rise2XNeg = rise2XNeg;
        myStruct(j).rise2VAll = rise2VAll;
        myStruct(j).rise2VPos = rise2VPos;
        myStruct(j).rise2VNeg = rise2VNeg;
        myStruct(j).FWHMall = FWHM_vecAll;
        myStruct(j).FWHMpos = FWHM_vecPos;
        myStruct(j).FWHMneg = FWHM_vecNeg;
        
        myStruct(j).Vresiduals_matAll = Vresiduals_matAll;
        myStruct(j).Vresiduals_matPos = Vresiduals_matPos;
        myStruct(j).Vresiduals_matNeg = Vresiduals_matNeg;
        myStruct(j).tauResiduals_mat = tauResiduals_mat;
        myStruct(j).Name = Name;
        
        myStruct(j).Events = relevantEvents;
        myStruct(j).residual_locs = residual_locs;
        
        myStruct(j).shifts_vec = shifts_cell;
        myStruct(j).RMS_vec = RMS_cell;
        
        myStruct(j).maxXi_loc_vec = maxXi_loc_vec;
        myStruct(j).maxYi_loc_vec = maxYi_loc_vec;
        myStruct(j).maxXi_vel_vec = maxXi_vel_vec;
        myStruct(j).maxYi_vel_vec = maxYi_vel_vec;
        
    catch
        disp([D{F(j)},' Failed to save']);
    end
end

save('DAandREsidualsStruct.mat','myStruct','-v7.3');