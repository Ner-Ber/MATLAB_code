function Movie_PlotRowFromHeatMap(RelevantRowOverTime,varargin)
% plot some frames out of a RowOverTime map. either use datatip and a
% single input to the function or u can enter the frame numbers (two input
% arguments).

if isempty(varargin{1})
    [~,y] = ginput;
    y = round(y);
else
    y = varargin{1};
end
row2plot = RelevantRowOverTime(y,:);
colors2Plot = MyVaryColor(length(y));
LGND = {};
figure; hold on;
for i = 1:length(y)
    plot(row2plot(i,:),'Color',colors2Plot(i,:));
    LGND{i} = ['row ',num2str(y(i))];
end
legend(LGND);
hold off;
end