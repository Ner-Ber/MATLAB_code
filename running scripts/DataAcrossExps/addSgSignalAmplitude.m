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
        UxxValues = zeros(3,N,length(relevantEvents));
        UxxMaxAmpAndLoc = zeros(2,N,length(relevantEvents));
        UyyValues = zeros(3,N,length(relevantEvents));
        UyyMaxAmpAndLoc = zeros(2,N,length(relevantEvents));
        UxyValues = zeros(3,N,length(relevantEvents));
        UxyMaxAmpAndLoc = zeros(2,N,length(relevantEvents));
        
        relvntRange = [-0.1 0.1];
        
        for iEv = 1:length(relevantEvents)
            try
                
                PhediStructCell{relevantEvents(iEv)}.PhediDataSG = PhediStructCell{relevantEvents(iEv)}.PhediData;
                PhediStructCell{relevantEvents(iEv)}.PhediData = PhediStructCell{relevantEvents(iEv)}.PhediDataSimpSmth;
                DataStruct = PhediStructCell{relevantEvents(iEv)};
                
                for Nsg = 1:15
                    %--- Uxx
                    XData_t = DataStruct.SgData.x_mins_x_tips(:,Nsg);
                    YData_t = 1e-3*(DataStruct.SgData.Uxx(:,Nsg)-mean(DataStruct.SgData.Uxx(1:100,Nsg)));
                    relvntLogic = XData_t>=min(relvntRange) & XData_t<=max(relvntRange);
                    XData = XData_t(relvntLogic);
                    YData = YData_t(relvntLogic);
                    
                    %- smooth ununiform data:
                    [x_out, ~, smth_out] = my_movVar(XData,YData,3e-3,0.5e-3);
                    %- find values a -0.05, 0.05, 0:
                    [~,I] = min(abs(bsxfun(@minus, x_out(:),[-0.05 0 0.05])));
                    ValuesUxx = smth_out(I);
                    [minValUxx,I] = min(smth_out);
                    minValLocUxx = x_out(I);
                    %- add data
                    UxxValues(:,Nsg,iEv) = ValuesUxx(:);
                    UxxMaxAmpAndLoc(:,Nsg,iEv) = [minValUxx; minValLocUxx];
                    
                    %--- Uyy
                    XData_t = DataStruct.SgData.x_mins_x_tips(:,Nsg);
                    YData_t = 1e-3*(DataStruct.SgData.Uyy(:,Nsg)-mean(DataStruct.SgData.Uyy(1:100,Nsg)));
                    relvntLogic = XData_t>=min(relvntRange) & XData_t<=max(relvntRange);
                    XData = XData_t(relvntLogic);
                    YData = YData_t(relvntLogic);
                    %- smooth ununiform data:
                    [x_out, ~, smth_out] = my_movVar(XData,YData,3e-3,0.5e-3);
                    %- find values a -0.05, 0.05, 0:
                    [~,I] = min(abs(bsxfun(@minus, x_out(:),[-0.05 0 0.05])));
                    ValuesUyy = smth_out(I);
                    [minValUyy,I] = max(smth_out);
                    minValLocUyy = x_out(I);
                    %- add data
                    UyyValues(:,Nsg,iEv) = ValuesUyy(:);
                    UyyMaxAmpAndLoc(:,Nsg,iEv) = [minValUyy; minValLocUyy];
                    
                    
                    %--- Uxy
                    zeroRegion = [-0.05 -0.065];
                    referenceSG_logical = (DataStruct.SgData.x_mins_x_tips(:,Nsg)<=zeroRegion(1)) & (DataStruct.SgData.x_mins_x_tips(:,Nsg)>=zeroRegion(2));
                    XData_t = DataStruct.SgData.x_mins_x_tips(:,Nsg);
                    YData_t = 1e-3*(DataStruct.SgData.Uxy(:,Nsg)-mean(DataStruct.SgData.Uxy(referenceSG_logical,Nsg)));
                    relvntLogic = XData_t>=min(relvntRange) & XData_t<=max(relvntRange);
                    XData = XData_t(relvntLogic);
                    YData = YData_t(relvntLogic);
                    %- smooth ununiform data:
                    [x_out, ~, smth_out] = my_movVar(XData,YData,3e-3,0.5e-3);
                    %- find values a -0.05, 0.05, 0:
                    [~,I] = min(abs(bsxfun(@minus, x_out(:),[-0.05 0 0.05])));
                    ValuesUxy = smth_out(I);
                    [minValUxy,I] = max(smth_out);
                    minValLocUxy = x_out(I);
                    %- add data
                    UxyValues(:,Nsg,iEv) = ValuesUxy(:);
                    UxyMaxAmpAndLoc(:,Nsg,iEv) = [minValUxy; minValLocUxy];
                    
                    
                end
            catch
                disp([D{F(j)},', ev=',num2str(relevantEvents(iEv)),' didn''t work']);
            end
        end
        Name = [PhediStructCell{relevantEvents(iEv)}.ExperimentData.ExpDate,' ',PhediStructCell{relevantEvents(iEv)}.ExperimentData.ExpHour];
        J = arrayfun(@(S) strcmpi(S.Name,Name), myStruct);
        
        myStruct(J).UxxValues = UxxValues;
        myStruct(J).UxxMaxAmpAndLoc = UxxMaxAmpAndLoc;
        myStruct(J).UyyValues = UyyValues;
        myStruct(J).UyyMaxAmpAndLoc = UyyMaxAmpAndLoc;
        myStruct(J).UxyValues = UxyValues;
        myStruct(J).UxyMaxAmpAndLoc = UxyMaxAmpAndLoc;
    catch
        disp([D{F(j)},' Failed to save']);
    end
end