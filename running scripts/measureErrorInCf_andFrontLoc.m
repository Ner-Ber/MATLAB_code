PatternEdges = [0.07 0.09];
relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
relevantEvents(relevantEvents==1) = [];
for j=1:length(relevantEvents)
    J = relevantEvents(j);
    PhediStructCell{J}.PhediDataSG = PhediStructCell{J}.PhediData;
    PhediStructCell{J}.PhediData = PhediStructCell{J}.PhediDataSimpSmth;
    
    BigPicRotStruct = PhediStructCell{J}.BigPicRotStruct;
    
    %% first plot pixel over time
    patternLogic = BigPicRotStruct.x>=min(PatternEdges) & BigPicRotStruct.x<=max(PatternEdges);
    pattern_pixles = BigPicRotStruct.DataMatNorm(:,patternLogic);
    
    ROT_time = BigPicRotStruct.t;
    
    x_plot = BigPicRotStruct.x(patternLogic);
    t_plot = BigPicRotStruct.frontTime_interp(patternLogic)/BigPicRotStruct.fps;
    
    A_front = nan(size(pattern_pixles,2),1);
    for i = 1:size(pattern_pixles,2)
        A_front(i) = interp1(ROT_time,pattern_pixles(:,i),t_plot(i),'linear');
    end
    
    V = [nan; diff(x_plot(:))./diff(t_plot(:))];
    
    Colors = MyVaryColor(size(pattern_pixles,2));
    
%     figure; hold on;
%     for i = 1:size(pattern_pixles,2)
%         plot(ROT_time,pattern_pixles(:,i),'Color',Colors(i,:));
%         plot(t_plot(i),A_front(i),'o','Color',Colors(i,:));
%     end
    
    %     figure; IDT_PlotRowOverTime(BigPicRotStruct);
    %     xlim([0.06 0.1]);
    %     figure; plot(V,A_front,'o');
    
    timeLogical = ROT_time>=min(t_plot) & ROT_time<=max(t_plot);
    pattern_pixles_space = BigPicRotStruct.DataMatNorm(timeLogical,:);
    cropped_time = BigPicRotStruct.t(timeLogical);
    A_front_space = nan(size(pattern_pixles_space,1),1);
    x_plot_space = nan(size(pattern_pixles_space,1),1);
    
    Colors = MyVaryColor(size(pattern_pixles_space,1));
%     figure; hold on;
    for i = 1:size(pattern_pixles_space,1)
        [~,I] = min(abs(cropped_time(i)-t_plot));
        x_plot_space(i) = x_plot(I);
        A_front_space(i) = interp1(BigPicRotStruct.x,pattern_pixles_space(i,:),x_plot_space(i),'linear');
        
%         plot(BigPicRotStruct.x,pattern_pixles_space(i,:),'Color',Colors(i,:));
%         plot(x_plot_space(i),A_front_space(i),'o','Color',Colors(i,:))

    end
%     ylim([0.5 1.2]);
    A_front_cell{j} = A_front;
    A_front_space_cell{j} = A_front_space;
    
end
