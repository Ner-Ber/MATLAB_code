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
        
        peakLoc_vecAll = zeros(length(relevantEvents),1);
        peakLoc_vecPos = zeros(length(relevantEvents),1);
        peakLoc_vecNeg = zeros(length(relevantEvents),1);
        for iEv = 1:length(relevantEvents)
            try
                
                PhediStructCell{relevantEvents(iEv)}.PhediDataSG = PhediStructCell{relevantEvents(iEv)}.PhediData;
                PhediStructCell{relevantEvents(iEv)}.PhediData = PhediStructCell{relevantEvents(iEv)}.PhediDataSimpSmth;
                
                BigPicRotStruct = PhediStructCell{relevantEvents(iEv)}.BigPicRotStruct;
                PhediData = PhediStructCell{relevantEvents(iEv)}.PhediData;
                slopeIncline = mean(PhediStructCell{relevantEvents(iEv)}.PhediData.slopeIncline,1);
                relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
                relevantEvents = relevantEvents(relevantEvents~=1);
                
                
                Sall = phedi_FWHMofMean(PhediStructCell{relevantEvents(iEv)},A_avgDist);
                VpeakXall = Sall.peakX;
                peakLoc_vecAll(iEv) = VpeakXall;
                
                
                posSlp = slopeIncline==+1;
                Spos = phedi_FWHMofMean(PhediStructCell{relevantEvents(iEv)},A_avgDist,posSlp);
                VpeakXpos = Spos.peakX;
                peakLoc_vecPos(iEv) = VpeakXpos;
                
                
                negSlp = slopeIncline==-1;
                Sneg = phedi_FWHMofMean(PhediStructCell{relevantEvents(iEv)},A_avgDist,negSlp);
                VpeakXneg = Sneg.peakX;
                peakLoc_vecNeg(iEv) = VpeakXneg;
                
            catch
            end
        end
        myStruct(j).peakLoc_vecAll = peakLoc_vecAll;
        myStruct(j).peakLoc_vecPos = peakLoc_vecPos;
        myStruct(j).peakLoc_vecNeg = peakLoc_vecNeg;
    catch
        disp([D{F(j)},' Failed to save']);
    end
end