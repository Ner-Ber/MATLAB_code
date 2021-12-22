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
        
        N = 15;
%         UyyResiduals_mat = zeros(length(residual_locs),N,length(relevantEvents));
%         UxxResiduals_mat = zeros(length(residual_locs),N,length(relevantEvents));
%         UxyResiduals_mat = zeros(length(residual_locs),N,length(relevantEvents));
        SyyResiduals_mat = zeros(length(residual_locs),N,length(relevantEvents));
        SxxResiduals_mat = zeros(length(residual_locs),N,length(relevantEvents));
        SxyResiduals_mat = zeros(length(residual_locs),N,length(relevantEvents));
%         Cf_at_sg = zeros(N,length(relevantEvents));
        
        for iEv = 1:length(relevantEvents)
            try
                
                PhediStructCell{relevantEvents(iEv)}.PhediDataSG = PhediStructCell{relevantEvents(iEv)}.PhediData;
                PhediStructCell{relevantEvents(iEv)}.PhediData = PhediStructCell{relevantEvents(iEv)}.PhediDataSimpSmth;
                
%                 %--- Uyy
%                 X_MIN_xTIP_sgOLD = PhediStructCell{relevantEvents(iEv)}.SgData.x_mins_x_tips;
%                 cropLogic = X_MIN_xTIP_sgOLD>=-0.1 & X_MIN_xTIP_sgOLD<=0.1;
%                 X_MIN_xTIP_sgOLD_cell = cellfun(@(A,B) A(B), num2cell(X_MIN_xTIP_sgOLD,1), num2cell(cropLogic,1),'UniformOutput',0);
%                 Uyy_old = bsxfun(@minus,PhediStructCell{relevantEvents(iEv)}.SgData.Uyy,mean(PhediStructCell{relevantEvents(iEv)}.SgData.Uyy(1:100,:)))*1e-3;
%                 Uyy_old_cropped = cellfun(@(A,B) A(B), num2cell(Uyy_old,1), num2cell(cropLogic,1),'UniformOutput',0);
%                 X_MIN_xTIP_sg = cellfun(@(A) max(A):-1e-4:min(A),X_MIN_xTIP_sgOLD_cell,'UniformOutput',0);
%                 Uyy = cellfun(@(A,B,C) interp1(A,B,C), X_MIN_xTIP_sgOLD_cell,Uyy_old_cropped,X_MIN_xTIP_sg,'UniformOutput',0);
%                 smthL = cellfun(@(A) A_avgDist/abs(diff(A([1,2]))), X_MIN_xTIP_sg,'UniformOutput',0);
%                 UyySmth = cellfun(@(A,B) smooth(A,B), Uyy,smthL,'UniformOutput',0);
%                 %- pad arrays with nans
%                 MV = max(cellfun(@length,UyySmth));
%                 UyySmthPaddRept = repmat(cell2mat(cellfun(@(A) padarray(A(:),[MV-length(A),0],nan,'post'),UyySmth,'UniformOutput',0)),1,1,length(residual_locs));
%                 X_MIN_xTIP_sgPadd = cell2mat(cellfun(@(A) padarray(A(:),[MV-length(A),0],nan,'post'),X_MIN_xTIP_sg,'UniformOutput',0));
%                 %- get minimums
%                 [~,I] = nanmin(abs(bsxfun(@minus,X_MIN_xTIP_sgPadd,repmat(permute(residual_locs(:),[3 2 1]),1,size(X_MIN_xTIP_sgPadd,2),1))));
%                 I2 = repmat(1:size(I,2),1,1,size(I,3));
%                 I3 = repmat(permute(1:size(I,3),[1 3 2]),1,size(I,2),1);
%                 II = sub2ind(size(UyySmthPaddRept),I(:),I2(:),I3(:));
%                 UyyResiduals = permute(reshape(UyySmthPaddRept(II),size(I)),[3 2 1]); % columns represent res locations, rows rep. sg num.
%                 UyyResiduals_mat(:,:,iEv) = UyyResiduals;
%                 
%                 %--- Uxx
%                 X_MIN_xTIP_sgOLD = PhediStructCell{relevantEvents(iEv)}.SgData.x_mins_x_tips;
%                 cropLogic = X_MIN_xTIP_sgOLD>=-0.1 & X_MIN_xTIP_sgOLD<=0.1;
%                 X_MIN_xTIP_sgOLD_cell = cellfun(@(A,B) A(B), num2cell(X_MIN_xTIP_sgOLD,1), num2cell(cropLogic,1),'UniformOutput',0);
%                 Uxx_old = bsxfun(@minus,PhediStructCell{relevantEvents(iEv)}.SgData.Uxx,mean(PhediStructCell{relevantEvents(iEv)}.SgData.Uxx(1:100,:)))*1e-3;
%                 Uxx_old_cropped = cellfun(@(A,B) A(B), num2cell(Uxx_old,1), num2cell(cropLogic,1),'UniformOutput',0);
%                 X_MIN_xTIP_sg = cellfun(@(A) max(A):-1e-4:min(A),X_MIN_xTIP_sgOLD_cell,'UniformOutput',0);
%                 Uxx = cellfun(@(A,B,C) interp1(A,B,C), X_MIN_xTIP_sgOLD_cell,Uxx_old_cropped,X_MIN_xTIP_sg,'UniformOutput',0);
%                 smthL = cellfun(@(A) A_avgDist/abs(diff(A([1,2]))), X_MIN_xTIP_sg,'UniformOutput',0);
%                 UxxSmth = cellfun(@(A,B) smooth(A,B), Uxx,smthL,'UniformOutput',0);
%                 %- pad arrays with nans
%                 MV = max(cellfun(@length,UxxSmth));
%                 UxxSmthPaddRept = repmat(cell2mat(cellfun(@(A) padarray(A(:),[MV-length(A),0],nan,'post'),UxxSmth,'UniformOutput',0)),1,1,length(residual_locs));
%                 X_MIN_xTIP_sgPadd = cell2mat(cellfun(@(A) padarray(A(:),[MV-length(A),0],nan,'post'),X_MIN_xTIP_sg,'UniformOutput',0));
%                 %- get minimums
%                 [~,I] = nanmin(abs(bsxfun(@minus,X_MIN_xTIP_sgPadd,repmat(permute(residual_locs(:),[3 2 1]),1,size(X_MIN_xTIP_sgPadd,2),1))));
%                 I2 = repmat(1:size(I,2),1,1,size(I,3));
%                 I3 = repmat(permute(1:size(I,3),[1 3 2]),1,size(I,2),1);
%                 II = sub2ind(size(UxxSmthPaddRept),I(:),I2(:),I3(:));
%                 UxxResiduals = permute(reshape(UxxSmthPaddRept(II),size(I)),[3 2 1]); % columns represent res locations, rows rep. sg num.
%                 UxxResiduals_mat(:,:,iEv) = UxxResiduals;
%                 
%                 
%                 %--- Uxy
%                 X_MIN_xTIP_sgOLD = PhediStructCell{relevantEvents(iEv)}.SgData.x_mins_x_tips;
%                 cropLogic = X_MIN_xTIP_sgOLD>=-0.1 & X_MIN_xTIP_sgOLD<=0.1;
%                 X_MIN_xTIP_sgOLD_cell = cellfun(@(A,B) A(B), num2cell(X_MIN_xTIP_sgOLD,1), num2cell(cropLogic,1),'UniformOutput',0);
%                 Uxy_old = bsxfun(@minus,PhediStructCell{relevantEvents(iEv)}.SgData.Uxy,mean(PhediStructCell{relevantEvents(iEv)}.SgData.Uxy(1:100,:)))*1e-3;
%                 Uxy_old_cropped = cellfun(@(A,B) A(B), num2cell(Uxy_old,1), num2cell(cropLogic,1),'UniformOutput',0);
%                 X_MIN_xTIP_sg = cellfun(@(A) max(A):-1e-4:min(A),X_MIN_xTIP_sgOLD_cell,'UniformOutput',0);
%                 Uxy = cellfun(@(A,B,C) interp1(A,B,C), X_MIN_xTIP_sgOLD_cell,Uxy_old_cropped,X_MIN_xTIP_sg,'UniformOutput',0);
%                 smthL = cellfun(@(A) A_avgDist/abs(diff(A([1,2]))), X_MIN_xTIP_sg,'UniformOutput',0);
%                 UxySmth = cellfun(@(A,B) smooth(A,B), Uxy,smthL,'UniformOutput',0);
%                 %- pad arrays with nans
%                 MV = max(cellfun(@length,UxySmth));
%                 UxySmthPaddRept = repmat(cell2mat(cellfun(@(A) padarray(A(:),[MV-length(A),0],nan,'post'),UxySmth,'UniformOutput',0)),1,1,length(residual_locs));
%                 X_MIN_xTIP_sgPadd = cell2mat(cellfun(@(A) padarray(A(:),[MV-length(A),0],nan,'post'),X_MIN_xTIP_sg,'UniformOutput',0));
%                 %- get minimums
%                 [~,I] = nanmin(abs(bsxfun(@minus,X_MIN_xTIP_sgPadd,repmat(permute(residual_locs(:),[3 2 1]),1,size(X_MIN_xTIP_sgPadd,2),1))));
%                 I2 = repmat(1:size(I,2),1,1,size(I,3));
%                 I3 = repmat(permute(1:size(I,3),[1 3 2]),1,size(I,2),1);
%                 II = sub2ind(size(UxySmthPaddRept),I(:),I2(:),I3(:));
%                 UxyResiduals = permute(reshape(UxySmthPaddRept(II),size(I)),[3 2 1]); % columns represent res locations, rows rep. sg num.
%                 UxyResiduals_mat(:,:,iEv) = UxyResiduals;
                
                
                %--- Syy
                X_MIN_xTIP_sgOLD = PhediStructCell{relevantEvents(iEv)}.SgData.x_mins_x_tips;
                cropLogic = X_MIN_xTIP_sgOLD>=-0.1 & X_MIN_xTIP_sgOLD<=0.1;
                X_MIN_xTIP_sgOLD_cell = cellfun(@(A,B) A(B), num2cell(X_MIN_xTIP_sgOLD,1), num2cell(cropLogic,1),'UniformOutput',0);
                Syy_old = PhediStructCell{relevantEvents(iEv)}.SgData.Syy;
                Syy_old_cropped = cellfun(@(A,B) A(B), num2cell(Syy_old,1), num2cell(cropLogic,1),'UniformOutput',0);
                X_MIN_xTIP_sg = cellfun(@(A) max(A):-1e-4:min(A),X_MIN_xTIP_sgOLD_cell,'UniformOutput',0);
                Syy = cellfun(@(A,B,C) interp1(A,B,C), X_MIN_xTIP_sgOLD_cell,Syy_old_cropped,X_MIN_xTIP_sg,'UniformOutput',0);
                smthL = cellfun(@(A) A_avgDist/abs(diff(A([1,2]))), X_MIN_xTIP_sg,'UniformOutput',0);
                SyySmth = cellfun(@(A,B) smooth(A,B), Syy,smthL,'UniformOutput',0);
                %- pad arrays with nans
                MV = max(cellfun(@length,SyySmth));
                SyySmthPaddRept = repmat(cell2mat(cellfun(@(A) padarray(A(:),[MV-length(A),0],nan,'post'),SyySmth,'UniformOutput',0)),1,1,length(residual_locs));
                X_MIN_xTIP_sgPadd = cell2mat(cellfun(@(A) padarray(A(:),[MV-length(A),0],nan,'post'),X_MIN_xTIP_sg,'UniformOutput',0));
                %- get minimums
                [~,I] = nanmin(abs(bsxfun(@minus,X_MIN_xTIP_sgPadd,repmat(permute(residual_locs(:),[3 2 1]),1,size(X_MIN_xTIP_sgPadd,2),1))));
                I2 = repmat(1:size(I,2),1,1,size(I,3));
                I3 = repmat(permute(1:size(I,3),[1 3 2]),1,size(I,2),1);
                II = sub2ind(size(SyySmthPaddRept),I(:),I2(:),I3(:));
                SyyResiduals = permute(reshape(SyySmthPaddRept(II),size(I)),[3 2 1]); % columns represent res locations, rows rep. sg num.
                SyyResiduals_mat(:,:,iEv) = SyyResiduals;
                
                %--- Sxx
                X_MIN_xTIP_sgOLD = PhediStructCell{relevantEvents(iEv)}.SgData.x_mins_x_tips;
                cropLogic = X_MIN_xTIP_sgOLD>=-0.1 & X_MIN_xTIP_sgOLD<=0.1;
                X_MIN_xTIP_sgOLD_cell = cellfun(@(A,B) A(B), num2cell(X_MIN_xTIP_sgOLD,1), num2cell(cropLogic,1),'UniformOutput',0);
                Sxx_old = PhediStructCell{relevantEvents(iEv)}.SgData.Sxx;
                Sxx_old_cropped = cellfun(@(A,B) A(B), num2cell(Sxx_old,1), num2cell(cropLogic,1),'UniformOutput',0);
                X_MIN_xTIP_sg = cellfun(@(A) max(A):-1e-4:min(A),X_MIN_xTIP_sgOLD_cell,'UniformOutput',0);
                Sxx = cellfun(@(A,B,C) interp1(A,B,C), X_MIN_xTIP_sgOLD_cell,Sxx_old_cropped,X_MIN_xTIP_sg,'UniformOutput',0);
                smthL = cellfun(@(A) A_avgDist/abs(diff(A([1,2]))), X_MIN_xTIP_sg,'UniformOutput',0);
                SxxSmth = cellfun(@(A,B) smooth(A,B), Sxx,smthL,'UniformOutput',0);
                %- pad arrays with nans
                MV = max(cellfun(@length,SxxSmth));
                SxxSmthPaddRept = repmat(cell2mat(cellfun(@(A) padarray(A(:),[MV-length(A),0],nan,'post'),SxxSmth,'UniformOutput',0)),1,1,length(residual_locs));
                X_MIN_xTIP_sgPadd = cell2mat(cellfun(@(A) padarray(A(:),[MV-length(A),0],nan,'post'),X_MIN_xTIP_sg,'UniformOutput',0));
                %- get minimums
                [~,I] = nanmin(abs(bsxfun(@minus,X_MIN_xTIP_sgPadd,repmat(permute(residual_locs(:),[3 2 1]),1,size(X_MIN_xTIP_sgPadd,2),1))));
                I2 = repmat(1:size(I,2),1,1,size(I,3));
                I3 = repmat(permute(1:size(I,3),[1 3 2]),1,size(I,2),1);
                II = sub2ind(size(SxxSmthPaddRept),I(:),I2(:),I3(:));
                SxxResiduals = permute(reshape(SxxSmthPaddRept(II),size(I)),[3 2 1]); % columns represent res locations, rows rep. sg num.
                SxxResiduals_mat(:,:,iEv) = SxxResiduals;
                
                
                %--- Sxy
                X_MIN_xTIP_sgOLD = PhediStructCell{relevantEvents(iEv)}.SgData.x_mins_x_tips;
                cropLogic = X_MIN_xTIP_sgOLD>=-0.1 & X_MIN_xTIP_sgOLD<=0.1;
                X_MIN_xTIP_sgOLD_cell = cellfun(@(A,B) A(B), num2cell(X_MIN_xTIP_sgOLD,1), num2cell(cropLogic,1),'UniformOutput',0);
                Sxy_old = PhediStructCell{relevantEvents(iEv)}.SgData.Sxy;
                Sxy_old_cropped = cellfun(@(A,B) A(B), num2cell(Sxy_old,1), num2cell(cropLogic,1),'UniformOutput',0);
                X_MIN_xTIP_sg = cellfun(@(A) max(A):-1e-4:min(A),X_MIN_xTIP_sgOLD_cell,'UniformOutput',0);
                Sxy = cellfun(@(A,B,C) interp1(A,B,C), X_MIN_xTIP_sgOLD_cell,Sxy_old_cropped,X_MIN_xTIP_sg,'UniformOutput',0);
                smthL = cellfun(@(A) A_avgDist/abs(diff(A([1,2]))), X_MIN_xTIP_sg,'UniformOutput',0);
                SxySmth = cellfun(@(A,B) smooth(A,B), Sxy,smthL,'UniformOutput',0);
                %- pad arrays with nans
                MV = max(cellfun(@length,SxySmth));
                SxySmthPaddRept = repmat(cell2mat(cellfun(@(A) padarray(A(:),[MV-length(A),0],nan,'post'),SxySmth,'UniformOutput',0)),1,1,length(residual_locs));
                X_MIN_xTIP_sgPadd = cell2mat(cellfun(@(A) padarray(A(:),[MV-length(A),0],nan,'post'),X_MIN_xTIP_sg,'UniformOutput',0));
                %- get minimums
                [~,I] = nanmin(abs(bsxfun(@minus,X_MIN_xTIP_sgPadd,repmat(permute(residual_locs(:),[3 2 1]),1,size(X_MIN_xTIP_sgPadd,2),1))));
                I2 = repmat(1:size(I,2),1,1,size(I,3));
                I3 = repmat(permute(1:size(I,3),[1 3 2]),1,size(I,2),1);
                II = sub2ind(size(SxySmthPaddRept),I(:),I2(:),I3(:));
                SxyResiduals = permute(reshape(SxySmthPaddRept(II),size(I)),[3 2 1]); % columns represent res locations, rows rep. sg num.
                SxyResiduals_mat(:,:,iEv) = SxyResiduals;


%                 %--- cf at the sg
%                 Cf_at_sg(:,iEv) = PhediStructCell{relevantEvents(iEv)}.SgData.v_sg(:);
                
                
            catch
                disp([D{F(j)},', ev=',num2str(relevantEvents(iEv)),' didn''t work']);
            end
        end
        Name = [PhediStructCell{relevantEvents(iEv)}.ExperimentData.ExpDate,' ',PhediStructCell{relevantEvents(iEv)}.ExperimentData.ExpHour];
        J = arrayfun(@(S) strcmpi(S.Name,Name), myStruct);
        
%         myStruct(J).UyyResiduals_mat = UyyResiduals_mat;
%         myStruct(J).UxxResiduals_mat = UxxResiduals_mat;
%         myStruct(J).UxyResiduals_mat = UxyResiduals_mat;
        myStruct(J).SyyResiduals_mat = SyyResiduals_mat;
        myStruct(J).SxxResiduals_mat = SxxResiduals_mat;
        myStruct(J).SxyResiduals_mat = SxyResiduals_mat;
%         myStruct(J).cfAtSg= Cf_at_sg;
    catch
        disp([D{F(j)},' Failed to save']);
    end
end