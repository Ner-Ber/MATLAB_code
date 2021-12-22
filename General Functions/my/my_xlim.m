function my_xlim(x_lim,fig_num_vec)
%The function changes the xlimmit for required figures, then chnges ylimmit
%to best fit the axes.
%Very usefull for changing axes to plotyy!!


if nargin<2
    fig_num_vec=gcf;
end


for j=1:length(fig_num_vec)
    ax_mesh=findobj(fig_num_vec(j),'Type','surface');
    if isempty(ax_mesh)

        %ax=get(fig_num_vec(j),'Children');
        ax=findobj(fig_num_vec(j),'Type','axes');
        set(ax,'Xlim',x_lim);
        set(ax,'XTickmode','auto');
        set(ax,'Ylimmode','auto');
        set(ax,'YTickmode','auto');
        
        if length(ax)>1 %special treatment for plotyy
            
            set(ax,'box','off')%'box' is looking bad with plotyy
            axx=cell2mat(get(ax,'Children'));
            set(axx,'Marker','.');
        end
        
    else
        ax_mesh=findobj(fig_num_vec(j),'Type','axes');
        ylim(ax_mesh(2),x_lim);
    end
end


