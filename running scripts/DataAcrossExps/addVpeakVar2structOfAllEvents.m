% D = dir_names
k = strfind(D,'.mat')
F = find(~cellfun(@isempty,k))
A_avgDist = 0.01;
for j=1:length(F)
    %% iterate on events
    try
        disp(['loading ',D{F(j)}]);
        load(D{F(j)});
        relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
        relevantEvents = relevantEvents(relevantEvents~=1);
        
        VarVPeak_vecAll = zeros(length(relevantEvents),1);
        for iEv = 1:length(relevantEvents)
            try
                
                PhediStructCell{relevantEvents(iEv)}.PhediDataSG = PhediStructCell{relevantEvents(iEv)}.PhediData;
                PhediStructCell{relevantEvents(iEv)}.PhediData = PhediStructCell{relevantEvents(iEv)}.PhediDataSimpSmth;
                
                relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
                relevantEvents = relevantEvents(relevantEvents~=1);
                
                
                Sall = phedi_FWHMofMean(PhediStructCell{relevantEvents(iEv)},A_avgDist);
                %-- get var at peak
                peakX = Sall.peakX;
                [~,I] = min(abs(Sall.AvgPhediStruct.X_mean-peakX));
                PeakVelVar = Sall.AvgPhediStruct.Vel_var_x(I);
                VarVPeak_vecAll(iEv) = PeakVelVar;
                
            catch
            end
        end
        myStruct(j).VarVpeak_vecAll = VarVPeak_vecAll;
    catch
        disp([D{F(j)},' Failed to save']);
    end
end