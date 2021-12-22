function my_plot_set_figure(fig_num)

if nargin==0
    fig_num=gcf; %gcf is handle for the current plot. if you want you can enter fig_num=1, i.e, apply on figure 1.
end

%----------Get handles
ax=findobj(fig_num,'Type','axes');
ax_mesh=findobj(fig_num,'Type','surface');
axx=get(ax,'Children');


%--set names
%ylabel(ax(end),'Left axis name','fontsize',30)
%ylabel(ax(1),'Right axis name','fontsize',30)

%---Plotyy: It's common that the two axes are not aligned. It worth trying
%some thing like that:

%1.set(ax(1),'ActivePositionProperty','position')
%2.set(ax(2),'ActivePositionProperty','position')
%3.set(ax(1),'position',get(ax(2),'position'))

%----change Font Size
%set(ax,'box','off')%'box' is looking bad with plotyy
for j=1:length(ax)
set(ax,'Fontsize',70)
%set(get(ax(j),'Title'),'Fontsize',40);
set(get(ax(j),'Xlabel'),'Fontsize',70);
set(get(ax(j),'Ylabel'),'Fontsize',70);
%set(ax(j),'Outerpos',[0 0 1 1]);
%when problems with the two axes on plotyy use - set(ax(1),'position',get(ax(2),'position'));
end

%------change Markers and lines
ax_mesh=findobj(fig_num,'Type','surface');
if isempty(ax_mesh)
    
    if iscell(axx)
        axx=cell2mat(axx);
    end
    %set(axx,'Marker','.');
    set(axx,'markerSize',40)% change marker Size %used 15
    %set(axx,'line','-'); % use 'none' fore no line
    %set(axx,'LineWidth',0.5)% %Line width
end

