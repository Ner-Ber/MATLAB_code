function my_plot_load_and_fix_fig(path_dir,name)

open([path_dir '\' name]);
%----------Get handles
ax=findobj(gcf,'Type','axes');
ax_mesh=findobj(gcf,'Type','surface');
axx=get(ax,'Children');


%----change Font Size
%set(ax,'box','off')%'box' is looking bad with plotyy
for j=1:length(ax)
set(ax,'Fontsize',50)
%set(get(ax(j),'Title'),'Fontsize',40);
set(get(ax(j),'Xlabel'),'Fontsize',50);
set(get(ax(j),'Ylabel'),'Fontsize',50);
%set(ax(j),'Outerpos',[0 0 1 1]);
end
pause
saveas(gcf,[path_dir '\' 'try.fig'],'fig')

