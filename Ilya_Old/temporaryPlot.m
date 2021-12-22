%% general
EventsVec = 2:length(analyzePhediCell);
for Event = EventsVec
    PhediLocation = analyzePhediCell{Event}.PhediLocation;
    colormap jet;
    FigColors = MyVaryColor(size(PhediLocation,2));
    
    %% plot only phedi
    
    figure; hold on;
    for i = 1:size(PhediLocation,2)
        if analyzePhediCell{Event}.slopeIncline(i)<0
            plot(analyzePhediCell{Event}.timeVec,PhediLocation(:,i),...
                'x-','Color',FigColors(i,:));
        end
    end
    legendLogical = analyzePhediCell{Event}.slopeIncline<0;
    LGND = cellfun(@num2str,num2cell(analyzePhediCell{Event}.measuredPhedisFromPlot(legendLogical)),'UniformOutput',0);
    xlabel('tims [s]');
    ylabel('location of phedi [m]');
    
    %% plot sg
    relevantSGs = [7,9];
    colormap default;
    FigColors = MyVaryColor(size(sg_data_cell{Event}.Uxx,2));
    yyaxis right
    for i = relevantSGs
        plot(sg_data_cell{Event}.t./1000,sg_data_cell{Event}.Uxx(:,i),...
            '.-','Color',FigColors(i,:));
    end
    ylabel('Uxx');
    LGND = [LGND;cellfun(@num2str,num2cell(sg_data_cell{Event}.x_sg(relevantSGs)'),'UniformOutput',0)];
    legend(LGND,'Location','eastoutside');
    
    Cf = sg_calc_velocity_from_sg(sg_data_cell{Event}, relevantSGs);
    title(['Cf = ',num2str(Cf)]);
    xlim([-1 1]*2e-4);
    hold off;
    
    %% plot speed from sg
    %--- plot speed of phedis:
    PhediLocation = analyzePhediCell{Event}.PhediLocation;
    colormap jet;
    FigColors = MyVaryColor(size(PhediLocation,2));
    figure; hold on;
    for i = 1:size(PhediLocation,2)
        if analyzePhediCell{Event}.slopeIncline(i)<0
            currentPhediVelocity = diff(PhediLocation(:,i))'./diff(analyzePhediCell{Event}.timeVec);
            currentTime = movmean(analyzePhediCell{Event}.timeVec,2,'Endpoints','discard');
            plot(currentTime,currentPhediVelocity,...
                'x-','Color',FigColors(i,:));
        end
    end
    legendLogical = analyzePhediCell{Event}.slopeIncline<0;
    LGND = cellfun(@num2str,num2cell(analyzePhediCell{Event}.measuredPhedisFromPlot(legendLogical)),'UniformOutput',0);
    xlabel('tims [s]');
    ylabel('velocity of phedi [m/s]');
    ylim([-2, 2]);
    
    %--- plot particle velocity from sg
    current_Cf = sg_data_cell{Event}.Cf;
    relevantSGs = [7,9];
    colormap default;
    FigColors = MyVaryColor(size(sg_data_cell{Event}.Uxx,2));
    
    yyaxis right
    for i = relevantSGs
        CurrentParticleVelocity = -current_Cf*(sg_data_cell{Event}.Uxx(:,i)-mean(sg_data_cell{Event}.Uxx(1:100,i)))*1e-3;
        plot(sg_data_cell{Event}.t./1000,CurrentParticleVelocity,...
            '.-','Color',FigColors(i,:));
    end
    ylabel('particle velocity (m/s)');
    LGND = [LGND;cellfun(@num2str,num2cell(sg_data_cell{Event}.x_sg(relevantSGs)'),'UniformOutput',0)];
    legend(LGND,'Location','eastoutside');
    
   
    title({['Event ',num2str(Event)],...
        ['Cf = ',num2str(current_Cf),' [m/s]'],...
        ['v_{phedi} ',num2str(analyzePhediCell{Event}.Phedi_velocity),' [m/s]']});
    xlim([-1 1]*2e-4);
    hold off;
end

