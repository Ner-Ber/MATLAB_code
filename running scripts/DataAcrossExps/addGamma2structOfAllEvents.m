% D = dir_names;
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
        Gamma_vec = zeros(length(relevantEvents),1);
        for iEv = 1:length(relevantEvents)
            try
                Gamma_vec(iEv) = PhediStructCell{relevantEvents(iEv)}.solAtSG.Gamma;
            catch
            end
        end
        myStruct(j).Gamma = Gamma_vec;
    catch
        disp([D{F(j)},' Failed to save']);
    end
end