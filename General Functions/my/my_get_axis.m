function fig=my_get_axis(fig_num,line_num)
%if using plotyy you may get the data from specific axes. Just enter
%instead of fig_num the axes handle.
%example:  ax=findobj(fig_num,'Type','axes'); -> fig=my_get_axis(ax(1))
%
%the function gets only visible lines on the plot unless "lin_num" is used.
%In this case the function gets the asked lines from thos which are visible

%-other examples
%a=get(gca,'Children'); c=get(a,'color'); set(a,'color',[ r g b]); set(a,'color',c{j})
if nargin<1
fig_num=gca;
else
    figure(fig_num);
end
ax=findobj(fig_num,'Type','surface');
%ax=findobj(fig_num,'Type','patch');


if isempty(ax)
    %----get axis from plot
    %axe=findobj(fig_num,'Type','line'); Instead this command,get first
    %axes handle and then the 'line' handle. needed if using plotyy
    
    %ax=findobj(fig_num,'Type','axes');%if fig_num is axes handle the command returns ax=fig_num
    ax=findobj(gca,'Type','axes');
    axe=findobj(ax,'Type','line');
    
    %-get only selected lines
        visbleVec=get(axe,'visible');
        visbleVec=strcmp( visbleVec,'on');
        axe=axe(visbleVec);
    if nargin>1
        axe=axe(length(axe)-line_num+1);
    end
    
    x=get(axe,'XData');   %cell array of the data
    y=get(axe,'YData');
    DisplayName=get(axe,'DisplayName');
    c=get(axe,'color'); %To change color use these: set(axe{j},'color',[r g b]); set(axe{j},'color',c{j})
    fig.title=get(get(ax,'title'),'string');
    fig.xLabel=get(get(ax,'xLabel'),'string');
    fig.yLabel=get(get(ax,'yLabel'),'string');
    fig.XLim=get(ax,'XLim');
    fig.YLim=get(ax,'YLim');
    fig.XTick=get(ax,'XTick');
    fig.YTick=get(ax,'YTick');
    
    if iscell(x)
        x=x(end:-1:1,1); %rearange the data -> left column is first plot
        y=y(end:-1:1,1);
        DisplayName=DisplayName(end:-1:1,1);
        c=c(end:-1:1,1);
        
        for j=1:length(x) %must be anicer way to do it
            x_size(j)=length(x{j});
        end
        %--- if all cell arryas are the same size convert to mat
        if diff(x_size)==0
            x=cell2mat(x)';
            y=cell2mat(y)';
            %--if the plots have same time axis return only one
            if diff(x,1,2)==0;
                x=x(:,1);
            end
        end
    else
        x=x';% Return column vector
        y=y';
    end
    fig.x=x;
    fig.y=y;
    fig.DisplayName=DisplayName;
    fig.color=c;
    %----get axis from mesh
else
    x=get(ax,'XData');
    t=get(ax,'YData');
    fig.z=get(ax,'ZData');
    fig.x=x(1,:);            %when geting data from mesh the axis are matrices, convert them back to vectors
    fig.t=t(:,1);
end



