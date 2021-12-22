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
        
        
        SyyResiduals_mat = zeros(length(residual_locs),length(relevantEvents));
        SxxResiduals_mat = zeros(length(residual_locs),length(relevantEvents));
        
        for iEv = 1:length(relevantEvents)
            try
                
                PhediStructCell{relevantEvents(iEv)}.PhediDataSG = PhediStructCell{relevantEvents(iEv)}.PhediData;
                PhediStructCell{relevantEvents(iEv)}.PhediData = PhediStructCell{relevantEvents(iEv)}.PhediDataSimpSmth;
                
                %--- Syy
                X_MIN_xTIP_sgOLD = PhediStructCell{relevantEvents(iEv)}.SgData.x_mins_x_tips(:,8);
                cropLogic = X_MIN_xTIP_sgOLD>=-0.1 & X_MIN_xTIP_sgOLD<=0.1;
                X_MIN_xTIP_sgOLD = X_MIN_xTIP_sgOLD(cropLogic);
                Syy_old = PhediStructCell{relevantEvents(iEv)}.SgData.Syy(:,8);
                Syy_old = Syy_old(cropLogic);
                X_MIN_xTIP_sg = max(X_MIN_xTIP_sgOLD):-1e-4:min(X_MIN_xTIP_sgOLD);
                Syy = interp1(X_MIN_xTIP_sgOLD,Syy_old,X_MIN_xTIP_sg);
                smthL = A_avgDist/abs(diff(X_MIN_xTIP_sg([1,2])));
                SyySmth = smooth(Syy,smthL);
                [~,I] = min(abs(bsxfun(@minus,X_MIN_xTIP_sg(:),residual_locs(:)')));
                SyyResiduals = SyySmth(I);
                SyyResiduals_mat(:,iEv) = SyyResiduals(:);
                
                %--- Sxx
                X_MIN_xTIP_sgOLD = PhediStructCell{relevantEvents(iEv)}.SgData.x_mins_x_tips(:,8);
                cropLogic = X_MIN_xTIP_sgOLD>=-0.1 & X_MIN_xTIP_sgOLD<=0.1;
                X_MIN_xTIP_sgOLD = X_MIN_xTIP_sgOLD(cropLogic);
                Sxx_old = PhediStructCell{relevantEvents(iEv)}.SgData.Sxx(:,8);
                Sxx_old = Sxx_old(cropLogic);
                X_MIN_xTIP_sg = max(X_MIN_xTIP_sgOLD):-1e-4:min(X_MIN_xTIP_sgOLD);
                Sxx = interp1(X_MIN_xTIP_sgOLD,Sxx_old,X_MIN_xTIP_sg);
                smthL = A_avgDist/abs(diff(X_MIN_xTIP_sg([1,2])));
                SxxSmth = smooth(Sxx,smthL);
                [~,I] = min(abs(bsxfun(@minus,X_MIN_xTIP_sg(:),residual_locs(:)')));
                SxxResiduals = SxxSmth(I);
                SxxResiduals_mat(:,iEv) = SxxResiduals(:);
                
                
            catch
            end
        end
        Name = [PhediStructCell{relevantEvents(iEv)}.ExperimentData.ExpDate,' ',PhediStructCell{relevantEvents(iEv)}.ExperimentData.ExpHour];
        J = arrayfun(@(S) strcmpi(S.Name,Name), myStruct);
        
        myStruct(J).SyyResiduals_mat = SyyResiduals_mat;
        myStruct(J).SxxResiduals_mat = SxxResiduals_mat;
    catch
        disp([D{F(j)},' Failed to save']);
    end
end