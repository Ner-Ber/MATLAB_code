function HeatMap = PlotRowOverTime(RowOverTime, fps, pixel_size)

MatSize = size(RowOverTime);
frameDuration = 1/fps;

x = linspace(0,MatSize(2)*pixel_size,MatSize(2));
y = linspace(0,MatSize(1)*frameDuration,MatSize(1));
HeatMap = imagesc('XData',x,'YData',y,'CData',RowOverTime);
xlim([0 x(end)]);  ylim([0 y(end)]);


end