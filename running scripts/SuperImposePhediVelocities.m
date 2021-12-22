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
        CfRegions = linspace(0,1255,100);
        CfVals = movmean(CfRegions,2,'Endpoints','discard');
        Colors = MyVaryColor(length(CfVals),flipud(My_colorMap));
        figure; ax=axes; hold on;
        maxCf = -inf;
        for i=1:length(relevantEvents)
            DataStruct = PhediStructCell{relevantEvents(i)};
            %-- load mean phedi
            if strcmpi(sideStr,'neg');
                phdIdx = mean(DataStruct.PhediData.slopeIncline,1)==-1;
            elseif strcmpi(sideStr,'pos');
                phdIdx = mean(DataStruct.PhediData.slopeIncline,1)==+1;
            else
                phdIdx = ones(1,size(DataStruct.PhediData.slopeIncline,2));
            end
            AvgPhediStruct = phedi_averagePhedi(DataStruct,phdIdx);
            x_Phedi = AvgPhediStruct.X_mean;
            vel_phedi = AvgPhediStruct.Vel_mean_x;
            
            %% plot
            [~,I] = min(abs(CfVals-abs(DataStruct.PhediData.Cf)));
            Name = [DataStruct.BigPicRotStruct.details,' Cf=',num2str(DataStruct.PhediData.Cf)];
            plot(x_Phedi,vel_phedi,'Color',Colors(I,:),'DisplayName',Name);
            xlabel('x_x_{tip} [m]');
            ylabel('phedi velocity [m/s]');
            
            maxCf = max(maxCf,max(vel_phedi));
        end
        plot([-0.05 0.05],[1 1]*0.4255,'--k');
        Title = [DataStruct.ExperimentData.ExpDate,' ',DataStruct.ExperimentData.ExpHour,' ',sideStr];
        title(Title);
        colormap(Colors);
        CB = colorbar;
        CB.Ticks = movmean(0:200:1255,2,'Endpoints','discard')/1255;
        CB.TickLabels = cellfun(@num2str,num2cell(round(CB.Ticks*1255)),'UniformOutput',0);
        axis([[-0.05 0.05]*1 [-0.1 maxCf*1.1]*1e-0]);
        saveCurrentFigure('C:\Users\NeriB\Google Drive\JAY_lab\THESIS\FIGURES\phediVelSuperimposed');
    catch
        disp(['unable ',D{k}]);
    end
end