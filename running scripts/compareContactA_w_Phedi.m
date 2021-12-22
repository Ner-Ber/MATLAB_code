%% compare contact area with phedi
% for k=1:length(D)
%     load(D{k});
    relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell)); relevantEvents = relevantEvents(relevantEvents~=1)
    for j=1:length(relevantEvents)
        PhediStructCell{relevantEvents(j)}.PhediDataSG = PhediStructCell{relevantEvents(j)}.PhediData;
        PhediStructCell{relevantEvents(j)}.PhediData = PhediStructCell{relevantEvents(j)}.PhediDataSimpSmth;
    end
    relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell)); relevantEvents = relevantEvents(relevantEvents~=1);
    CfRegions = 0:200:1200;
    CfVals = movmean(CfRegions,2,'Endpoints','discard');
    Colors = MyVaryColor(length(CfVals),parula);
    figure; ax=axes; hold on;
    for i=1:length(relevantEvents)
        DataStruct = PhediStructCell{relevantEvents(i)};
        %-- load mean phedi
        AvgPhediStruct = phedi_averagePhedi(DataStruct);
        x_Phedi = AvgPhediStruct.X_mean;
        vel_phedi = AvgPhediStruct.Vel_mean_x;
        %     vel_phedi = AvgPhediStruct.Loc_mean_x;
        
        %-- contact area
        thisLoc = DataStruct.PhediData.PhotoLocation;
        x_Contact = DataStruct.BigPicRotStruct.x - thisLoc;
        frontTime = DataStruct.BigPicRotStruct.frontTime_interp./DataStruct.BigPicRotStruct.fps;
        [~,I1] = min(abs(x_Contact));
        t_minTtip_4tmprl = DataStruct.BigPicRotStruct.t-frontTime(I1);
        [~,I2] = min(abs(t_minTtip_4tmprl));
        A_spatial = DataStruct.BigPicRotStruct.DataMatNorm(I2,:);
        
        %% create common axis
        X = [x_Contact(:);x_Phedi(:)];
        %--- interpolate for common spatial steps
        x_unique = unique(X);
        %-- exclude nans:
        DD = diff(x_Phedi);
        %     X_vec = min(x_unique(abs(x_unique)~=inf)):(abs(min(min(DD(abs(DD)~=inf))))/4):max(x_unique(abs(x_unique)~=inf));
        %     X_vec = X_vec(abs(X_vec)<=0.05);
        X_vec = linspace(min(x_unique(abs(x_unique)~=inf)),max(x_unique(abs(x_unique)~=inf)),5e4);
        X_vec = X_vec(abs(X_vec)<=0.05);
        %% interpolate on common axis
        vel_phedi_interp = interp1(x_Phedi,vel_phedi,X_vec);
        A_phedi_interp = interp1(x_Contact,A_spatial,X_vec);
        
        %% normalize velocity
        %     S = phedi_FWHMofMean(DataStruct);
        %     vel_phedi_interp = vel_phedi_interp./S.peakY;
        
        %% plot
        [~,I] = min(abs(CfVals-DataStruct.PhediData.Cf));
        Name = [DataStruct.BigPicRotStruct.details,' Cf=',num2str(DataStruct.PhediData.Cf)];
        plot3(X_vec,vel_phedi_interp,A_phedi_interp,'Color',Colors(I,:),'DisplayName',Name);
        xlabel('x_x_{tip} [m]');
        ylabel('phedi velocity [m/s]');
        zlabel('contact area [%]');
    end
    title([DataStruct.ExperimentData.ExpDate,' ',DataStruct.ExperimentData.ExpHour]);
    ax.View = [90 0];
    colormap(Colors);
    CB = colorbar;
    CB.Ticks = CfVals/1200;
    CB.TickLabels = cellfun(@num2str,num2cell(CfVals),'UniformOutput',0);
    saveCurrentFigure('C:\Users\NeriB\Google Drive\JAY_lab\THESIS\FIGURES\contactVSvelocity');
% end