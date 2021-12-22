function copyPlotStyle(originFigName, targetFigName)
% copyPlotStyle(originFig, targetFig)
%
% copy the style of each plot in the figure to a target figure. this will
% work well in the case where both figures have identical nu,ber of plots.
%
% in the case of entering only a origin fig, the fucntion will act on the
% current figure open

originFigHandle = originFigName;
GCA = get(get(originFigHandle,'Children'),'Children');
if nargin<2
    thisFig = get(gca,'Children');
else
    targetFigHandle = targetFigName;
    thisFig = get(get(targetFigHandle,'Children'),'Children');
end

for i=length(thisFig):-1:1
    thisFig(i).Color = GCA(i).Color;
    thisFig(i).Marker = GCA(i).Marker;
    thisFig(i).DisplayName = GCA(i).DisplayName;
end

end