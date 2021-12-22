function my_plot_change_axis_units(axis_name,ax_num,tick,tick_label,fig_num)
%if nargin==2 set relavant axes to auto


if nargin<5
    fig_num=gcf; %gcf is handle for the current plot. if you want you can enter fig_num=1, i.e, apply on figure 1.
end

%----------Get handles
ax=findobj(fig_num,'Type','axes');
ax_mesh=findobj(fig_num,'Type','surface');
axx=get(ax,'Children');

%----change x axis
%set(ax,'XticklabelMode','auto'),set(ax,'XtickMode','auto') return to original axis. 
% notice! xlim (ylim) is using 'Xtick' as refference position


set(ax(ax_num),[axis_name 'tickMode'],'auto')
set(ax(ax_num),[axis_name 'ticklabelMode'],'auto')

if nargin>2
set(ax(ax_num),[axis_name 'tick'],tick);
set(ax(ax_num),[axis_name 'ticklabel'],tick_label);
end


