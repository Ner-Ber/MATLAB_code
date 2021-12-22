%% superimpsoe phedivelocities
for k=1:length(D)
    load(D{k});
    try
        relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell)); relevantEvents = relevantEvents(relevantEvents~=1)
        for j=1:length(relevantEvents)
            PhediStructCell{relevantEvents(j)}.PhediDataSG = PhediStructCell{relevantEvents(j)}.PhediData;
            PhediStructCell{relevantEvents(j)}.PhediData = PhediStructCell{relevantEvents(j)}.PhediDataSimpSmth;
        end
        relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell)); relevantEvents = relevantEvents(relevantEvents~=1);
        GamRegions = linspace(0,10,10);
        GamaVals = movmean(GamRegions,2,'Endpoints','discard');
        Colors = MyVaryColor(length(GamaVals),flipud(My_colorMap));
        % figure; ax=axes; hold on;
        maxVpeak = -inf;
        for i=1:length(relevantEvents)
            DataStruct = PhediStructCell{relevantEvents(i)};
            %-- load mean phedi
            S = phedi_FWHMofMean(DataStruct);
            peakY = S.peakY;
            G = DataStruct.solAtSG.Gamma;
            Cf = DataStruct.PhediData.Cf;
            %% plot
            [~,I] = min(abs(GamaVals-G));
            Name = [DataStruct.BigPicRotStruct.details,' Cf=',num2str(DataStruct.PhediData.Cf),' \Gamma=',num2str(G)];
            plot(Cf,peakY,'o','Color',Colors(I,:),'DisplayName',Name);
            xlabel('C_{f} [m/s]');
            ylabel('v_{peak} [m/s]');
            
            maxVpeak = max(maxVpeak,max(peakY));
            drawnow;
        end
        title([DataStruct.ExperimentData.ExpDate,' ',DataStruct.ExperimentData.ExpHour]);
        colormap(Colors);
        CB = colorbar;
        CB.Ticks = GamaVals/10;
        CB.TickLabels = cellfun(@num2str,num2cell(round(GamaVals)),'UniformOutput',0);
        axis auto;
        % saveCurrentFigure('C:\Users\NeriB\Google Drive\JAY_lab\THESIS\FIGURES\phediVelSuperimposed');
    catch
        disp(['unable ',D{k}]);
    end
end