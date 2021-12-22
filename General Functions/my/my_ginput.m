function[x,y,line_list]=my_ginput(LineNum)
%Enter LineNum for only one graph. Other wise each graph will be opened
%seperatly . Note that each graph you take same amount of points.

if (nargin==0)
    LineNum='all';
end

plotbrowser('on');
ax=findobj(gca,'type','line');
set(ax,'Visible','off');
k=1;

for j=1:length(ax)
    
    if (strcmp(LineNum,'all')||(LineNum==length(ax)+1-j))
        
        set(ax(length(ax)+1-j),'Visible','on');
        [xx,yy]=ginput;
        set(ax(length(ax)+1-j),'Visible','off');
        if ~isempty(xx)
            x(:,k)=xx;
            y(:,k)=yy;
            line_list(k)=j;
            k=k+1;
        end
        
    end
end
set(ax,'Visible','on');