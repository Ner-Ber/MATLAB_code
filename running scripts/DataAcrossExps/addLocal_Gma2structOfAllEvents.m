load('Other\DAandREsidualsStruct.mat');
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
        
        LocalGammaVec = cell(length(relevantEvents),1);
        stretch_cell = cell(length(relevantEvents),1);
        RMS_strtch_cell = cell(length(relevantEvents),1);
        for iEv = 1:length(relevantEvents)
            try
                
                PhediStructCell{relevantEvents(iEv)}.PhediDataSG = PhediStructCell{relevantEvents(iEv)}.PhediData;
                PhediStructCell{relevantEvents(iEv)}.PhediData = PhediStructCell{relevantEvents(iEv)}.PhediDataSimpSmth;
                
                relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
                relevantEvents = relevantEvents(relevantEvents~=1);
                
                
                
                [~, RMS_strtch_vec, stretch_vec] = phedi_UxBestFitVerticalStretch(PhediStructCell{relevantEvents(iEv)});
                
                %-- local gamma vec
                fixedK = PhediStructCell{relevantEvents(iEv)}.solAtSG.K*stretch_vec;
                v = PhediStructCell{relevantEvents(iEv)}.PhediData.Cf;
                k=Cs/Cd;%Broberg p.330
                alpha_d=(1-(v./Cd).^2).^0.5;
                alpha_s=(1-(v./Cs).^2).^0.5;
                RlyghFunc=4*alpha_d.*alpha_s-(1+alpha_s.^2).^2;
                %----- Following Broberg p.334,336
                A=2*(1-k^2)*alpha_s.*v.^2./RlyghFunc/Cs^2;
                localG = (fixedK.^2)./(mu*4*(1-k^2)./A);
                LocalGammaVec{iEv} = localG;
                
                %                 [residualsCell, RMS_vec, shifts_vec] = phedi_UxBestFitVerticalShift(PhediStructCell{relevantEvents(iEv)});
                stretch_cell{iEv} = stretch_vec;
                RMS_strtch_cell{iEv} = RMS_strtch_vec;
                
            catch
                dusp('');
            end
        end
        
        Name = [PhediStructCell{relevantEvents(iEv)}.ExperimentData.ExpDate,' ',PhediStructCell{relevantEvents(iEv)}.ExperimentData.ExpHour];
        J = arrayfun(@(S) strcmpi(S.Name,Name), myStruct);
        
        myStruct(J).localGvec_vec = LocalGammaVec;
        myStruct(J).stretch_vec = stretch_cell;
        myStruct(J).RMS_strtch_vec = RMS_strtch_cell;
        
    catch
        disp([D{F(j)},' Failed to save']);
    end
end