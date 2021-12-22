% figure; hold on;
% colormap jet;
% Current_frames = 950:1200;
% FigColors = MyVaryColor(length(Current_frames));
% my_x_axis = 1:size(RowOverTimeCell{1},2);
% for i = 1:length(Current_frames)
%     plot(my_x_axis-PhediLocation(ismember(frameCount,Current_frames(i))),...
%         RowOverTimeCell{1}(Current_frames(i),:),...
%         'Color',FigColors(i,:));
%     
% end

%% --
relevanrRowOverTime = RowOverTimeCell{9};
%% calculate phedis
tic;
[frameCount, PhediLocation, slopeIncline] = Movie_calculate_continues_shift(relevanrRowOverTime, 850, 1200, measuredPhedisFromPlot);
toc;

%% plot heat map and phedis dynamics

FigColors = MyVaryColor(size(PhediLocation,2));
figure; hold on;
colormap default;
imagesc(relevanrRowOverTime);
colormap jet
for i = 1:size(PhediLocation,2)
    if slopeIncline(i)>0
        plot(PhediLocation(:,i)+measuredPhedisFromPlot(i),frameCount,...
            '*-','Color',FigColors(i,:));
    elseif slopeIncline(i)<0
        plot(PhediLocation(:,i)+measuredPhedisFromPlot(i),frameCount,...
            '.-','Color',FigColors(i,:));
    else
        plot(PhediLocation(:,i)+measuredPhedisFromPlot(i),frameCount,...
            'r--');
    end
end
LGND = cellfun(@num2str,num2cell(measuredPhedisFromPlot),'UniformOutput',0);
LGND = cat(2,{'row over time'},LGND);
legend(LGND);

%% plot only phedi
figure; hold on;
colormap jet
for i = 1:size(PhediLocation,2)
    if slopeIncline(i)>0
        plot(frameCount,PhediLocation(:,i),...
            '*-','Color',FigColors(i,:));
    elseif slopeIncline(i)<0
        plot(frameCount,PhediLocation(:,i),...
            '.-','Color',FigColors(i,:));
        else
        plot(frameCount,PhediLocation(:,i),...
            'r--');
    end
end
LGND = cellfun(@num2str,num2cell(measuredPhedisFromPlot),'UniformOutput',0);
legend(LGND);

%% plot crossection
Movie_PlotRowFromHeatMap(relevanrRowOverTime,850:1200);