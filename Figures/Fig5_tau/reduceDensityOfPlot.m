WantedDens = 0.01/150;
g = get(gca,'children');
f = get(gca);
% n = find(strcmp(arrayfun(@(A) A.DisplayName, g,'UniformOutput',0),'mean phedi'));
n = find(strcmp(arrayfun(@(A) A.DisplayName, g,'UniformOutput',0),'pix over time'));
XData_trimmed = g(n).XData(g(n).XData>=min(f.XAxis.Limits) & g(n).XData<=max(f.XAxis.Limits));
YData_trimmed = g(n).YData(g(n).XData>=min(f.XAxis.Limits) & g(n).XData<=max(f.XAxis.Limits));
D = min(abs(diff(XData_trimmed)));
L = length(XData_trimmed);
XData_reduced = XData_trimmed(1:round(WantedDens/D):L);
YData_reduced = YData_trimmed(1:round(WantedDens/D):L);

h = plot(XData_reduced,YData_reduced,'-','Color',g(n).Color,'LineWidth',g(n).LineWidth);
uistack(h,'bottom');