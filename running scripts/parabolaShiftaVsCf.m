%% parabola shifta vs Cf
%
% for k=1:length(D)
%     load(D{k});
%     try
relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell)); relevantEvents = relevantEvents(relevantEvents~=1)
for j=1:length(relevantEvents)
    PhediStructCell{relevantEvents(j)}.PhediDataSG = PhediStructCell{relevantEvents(j)}.PhediData;
    PhediStructCell{relevantEvents(j)}.PhediData = PhediStructCell{relevantEvents(j)}.PhediDataSimpSmth;
end
% figure; ax=axes; hold on;
maxVpeak = -inf;
for i=1:length(relevantEvents)
    try
        DataStruct = PhediStructCell{relevantEvents(i)};
        [residualsCell, RMS_vec, shifts_vec] = phedi_UxBestFitVerticalShift(DataStruct);
        shiftsVar = var(shifts_vec);
        G = DataStruct.solAtSG.Gamma;
        Cf = DataStruct.PhediData.Cf;
        %% plot
        Name = [DataStruct.BigPicRotStruct.details,' Cf=',num2str(DataStruct.PhediData.Cf)];
        plot(Cf,shiftsVar,'o','DisplayName',Name);
        xlabel('C_{f} [m/s]');
        ylabel('shift var [m^2]');
        
        maxVpeak = max(maxVpeak,max(shiftsVar));
        drawnow;
    catch
    end
end
% title([DataStruct.ExperimentData.ExpDate,' ',DataStruct.ExperimentData.ExpHour]);
axis auto;
% saveCurrentFigure('C:\Users\NeriB\Google Drive\JAY_lab\THESIS\FIGURES\phediVelSuperimposed');
%     catch
%         disp(['unable ',D{k}]);
%     end
% end