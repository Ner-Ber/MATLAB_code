D = dir_names
k1 = strfind(D,'.mat')
k2 = strfind(D,'allPhedi')
F = find(~cellfun(@isempty,k1) & ~cellfun(@isempty,k2))
A_avgDist = 0.01;
residual_locs = -0.1:0.01:-0.01;
for j=1:length(F)
    %% iterate on events
    try
        disp(['loading ',D{F(j)}]);
        load(D{F(j)});
        relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
        relevantEvents = relevantEvents(relevantEvents~=1);
        
        
        UyyResiduals_mat = zeros(length(residual_locs),length(relevantEvents));
        UxxResiduals_mat = zeros(length(residual_locs),length(relevantEvents));
        UxyResiduals_mat = zeros(length(residual_locs),length(relevantEvents));
        
        for iEv = 1:length(relevantEvents)
            try
                
                PhediStructCell{relevantEvents(iEv)}.PhediDataSG = PhediStructCell{relevantEvents(iEv)}.PhediData;
                PhediStructCell{relevantEvents(iEv)}.PhediData = PhediStructCell{relevantEvents(iEv)}.PhediDataSimpSmth;
                
                %--- Uyy
                X_MIN_xTIP_sgOLD = PhediStructCell{relevantEvents(iEv)}.SgData.x_mins_x_tips(:,8);
                cropLogic = X_MIN_xTIP_sgOLD>=-0.1 & X_MIN_xTIP_sgOLD<=0.1;
                X_MIN_xTIP_sgOLD = X_MIN_xTIP_sgOLD(cropLogic);
                Uyy_old = (PhediStructCell{relevantEvents(iEv)}.SgData.Uyy(:,8)-mean(PhediStructCell{relevantEvents(iEv)}.SgData.Uyy(1:100,8)))*1e-3;
                Uyy_old = Uyy_old(cropLogic);
                X_MIN_xTIP_sg = max(X_MIN_xTIP_sgOLD):-1e-4:min(X_MIN_xTIP_sgOLD);
                Uyy = interp1(X_MIN_xTIP_sgOLD,Uyy_old,X_MIN_xTIP_sg);
                smthL = A_avgDist/abs(diff(X_MIN_xTIP_sg([1,2])));
                UyySmth = smooth(Uyy,smthL);
                [~,I] = min(abs(bsxfun(@minus,X_MIN_xTIP_sg(:),residual_locs(:)')));
                UyyResiduals = UyySmth(I);
                UyyResiduals_mat(:,iEv) = UyyResiduals(:);
                
                %--- Uxx
                X_MIN_xTIP_sgOLD = PhediStructCell{relevantEvents(iEv)}.SgData.x_mins_x_tips(:,8);
                cropLogic = X_MIN_xTIP_sgOLD>=-0.1 & X_MIN_xTIP_sgOLD<=0.1;
                X_MIN_xTIP_sgOLD = X_MIN_xTIP_sgOLD(cropLogic);
                Uxx_old = (PhediStructCell{relevantEvents(iEv)}.SgData.Uxx(:,8)-mean(PhediStructCell{relevantEvents(iEv)}.SgData.Uxx(1:100,8)))*1e-3;
                Uxx_old = Uxx_old(cropLogic);
                X_MIN_xTIP_sg = max(X_MIN_xTIP_sgOLD):-1e-4:min(X_MIN_xTIP_sgOLD);
                Uxx = interp1(X_MIN_xTIP_sgOLD,Uxx_old,X_MIN_xTIP_sg);
                smthL = A_avgDist/abs(diff(X_MIN_xTIP_sg([1,2])));
                UxxSmth = smooth(Uxx,smthL);
                [~,I] = min(abs(bsxfun(@minus,X_MIN_xTIP_sg(:),residual_locs(:)')));
                UxxResiduals = UxxSmth(I);
                UxxResiduals_mat(:,iEv) = UxxResiduals(:);
                
                %--- Uxy
                X_MIN_xTIP_sgOLD = PhediStructCell{relevantEvents(iEv)}.SgData.x_mins_x_tips(:,8);
                cropLogic = X_MIN_xTIP_sgOLD>=-0.1 & X_MIN_xTIP_sgOLD<=0.1;
                X_MIN_xTIP_sgOLD = X_MIN_xTIP_sgOLD(cropLogic);
                Uxy_old = (PhediStructCell{relevantEvents(iEv)}.SgData.Uxy(:,8)-mean(PhediStructCell{relevantEvents(iEv)}.SgData.Uxy(1:100,8)))*1e-3;
                Uxy_old = Uxy_old(cropLogic);
                X_MIN_xTIP_sg = max(X_MIN_xTIP_sgOLD):-1e-4:min(X_MIN_xTIP_sgOLD);
                Uxy = interp1(X_MIN_xTIP_sgOLD,Uxy_old,X_MIN_xTIP_sg);
                smthL = A_avgDist/abs(diff(X_MIN_xTIP_sg([1,2])));
                UxySmth = smooth(Uxy,smthL);
                [~,I] = min(abs(bsxfun(@minus,X_MIN_xTIP_sg(:),residual_locs(:)')));
                UxyResiduals = UxySmth(I);
                UxyResiduals_mat(:,iEv) = UxyResiduals(:);
                
                
            catch
            end
        end
        Name = [PhediStructCell{relevantEvents(iEv)}.ExperimentData.ExpDate,' ',PhediStructCell{relevantEvents(iEv)}.ExperimentData.ExpHour];
        J = arrayfun(@(S) strcmpi(S.Name,Name), myStruct);
        
        myStruct(J).UyyResiduals_mat = UyyResiduals_mat;
        myStruct(J).UxxResiduals_mat = UxxResiduals_mat;
        myStruct(J).UxyResiduals_mat = UxyResiduals_mat;
    catch
        disp([D{F(j)},' Failed to save']);
    end
end