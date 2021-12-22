function fig=my_legend(legendV)

%-other examples
%a=get(gca,'Children'); c=get(a,'color'); set(a,'color',[ r g b]); set(a,'color',c{j})

%----get axis from plot
%axe=findobj(fig_num,'Type','line'); Instead this command,get first
%axes handle and then the 'line' handle. needed if using plotyy

fig_num=gcf;

ax=findobj(fig_num,'Type','axes');%if fig_num is axes handle the command returns ax=fig_num
axe=findobj(ax,'Type','line');

%-get only selected lines
visbleVec=get(axe,'visible');
visbleVec=strcmp( visbleVec,'on');
axe=axe(visbleVec);

for j=1:length(axe)
    set(axe(length(axe)-j+1),'DisplayName',legendV(j,:));
end





