D = dir_names
k1 = strfind(D,'.mat')
k2 = strfind(D,'allPhedi')
F = find(~cellfun(@isempty,k1) & ~cellfun(@isempty,k2))
A_avgDist = 0.01;
for j=1:length(F)
    %% iterate on events
    try
        disp(['loading ',D{F(j)}]);
        load(D{F(j)});
        relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
        relevantEvents = relevantEvents(relevantEvents~=1);
        
        shifts_cell = cell(length(relevantEvents),1);
        RMS_cell = cell(length(relevantEvents),1);
        for iEv = 1:length(relevantEvents)
            try
                
                PhediStructCell{relevantEvents(iEv)}.PhediDataSG = PhediStructCell{relevantEvents(iEv)}.PhediData;
                PhediStructCell{relevantEvents(iEv)}.PhediData = PhediStructCell{relevantEvents(iEv)}.PhediDataSimpSmth;
                
                relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
                relevantEvents = relevantEvents(relevantEvents~=1);
                
                
                
                [residualsCell, RMS_vec, shifts_vec] = phedi_UxBestFitVerticalShift(PhediStructCell{relevantEvents(iEv)});
                shifts_cell{iEv} = shifts_vec;
                RMS_cell{iEv} = RMS_vec;
                
            catch
            end
        end
        myStruct(j).shifts_vec = shifts_cell;
        myStruct(j).RMS_vec = RMS_cell;
        myStruct(j).G = PhediStructCell{relevantEvents(iEv)}.solAtSG.Gamma;
    catch
        disp([D{F(j)},' Failed to save']);
    end
end