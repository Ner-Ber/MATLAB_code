% load('Other\DAandREsidualsStruct.mat');
load('DAandREsidualsStruct.mat');
D = dir_names
k1 = strfind(D,'.mat')
k2 = strfind(D,'allPhedi')
F = find(~cellfun(@isempty,k1) & ~cellfun(@isempty,k2))
A_avgDist = 0.01;
[Cd, Cs, Cr, nu, ro, E, mu, Gamma, PlaneStrain, tau_p, Xc0]=CrackSolutionMaterialProperties;
for j=1:length(F)
    %% iterate on events
    try
        disp(['loading ',D{F(j)}]);
        load(D{F(j)});
        relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
        relevantEvents = relevantEvents(relevantEvents~=1);
        GammaOfAllSgs = nan(length(relevantEvents),15);
        
        for iEv = 1:length(relevantEvents)
            try
                
                PhediStructCell{relevantEvents(iEv)}.PhediDataSG = PhediStructCell{relevantEvents(iEv)}.PhediData;
                PhediStructCell{relevantEvents(iEv)}.PhediData = PhediStructCell{relevantEvents(iEv)}.PhediDataSimpSmth;
                
                relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
                relevantEvents = relevantEvents(relevantEvents~=1);
                
                [~,GammaVec] = gamma_calcGammaAlongInterface(PhediStructCell{relevantEvents(iEv)});
                GammaOfAllSgs(iEv,:) = GammaVec(:)';
                
            catch
                disp('');
            end
        end
        
        Name = [PhediStructCell{relevantEvents(iEv)}.ExperimentData.ExpDate,' ',PhediStructCell{relevantEvents(iEv)}.ExperimentData.ExpHour];
        J = arrayfun(@(S) strcmpi(S.Name,Name), myStruct);
        
        myStruct(J).GammaOfAllSgs = GammaOfAllSgs;
        
    catch
        disp([D{F(j)},' Failed to save']);
    end
end