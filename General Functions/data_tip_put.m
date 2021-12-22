function data_tip_put(dtip,col_list,line_list)
%The function get structur with x,y fields and plots the datatip
%col_list - the columns that sholud be plotted . example: you have event
%1,2,3 and want to plot only events 1,2 enter -col_list=[1,2].
%col_list= 'all' -plot all the dtip data
%line_list - the 'lines' (curves) on the figure that should be ploted on.
%for example if figure contains event 1,2,3, and dtip struct has only event
%1,3 you should skip curve 2.
%if col_list specified and line_list not, function assignes  line_list=col_list;


if strcmp(col_list,'all')
    col_list=1:length(dtip.x(1,:));
end

if nargin<3
if ~isfield(dtip,'line_list')
    line_list=col_list;
else
    line_list=dtip.line_list;
    line_list=line_list(col_list);
end
end

%---get the figure handles (line handle)
line_h=findobj(gcf,'Type','line');%all the graph handles
line_h=line_h(end:-1:1);

%It's recommended to create all the datatip first, then to set position.
%otherwise, some thing goes wrong

%---creat random DataTips
numrow=length(dtip.x(:,1));
for j=1:length(col_list)
    datatip_num=min([find(isnan(dtip.x(:,col_list(j))),1,'first')-1 numrow]); %Or number of rows or when NAN
    for k=1:datatip_num
        cursor_h(k,j)= createDatatip(datacursormode,line_h(line_list(j)));
    end
end

drawnow

%---set position to datatip
for j=1:length(col_list)
    datatip_num=min([find(isnan(dtip.x(:,col_list(j))),1,'first')-1 numrow]);
    for k=1:datatip_num
        set(cursor_h(k,j),'Position',[dtip.x(k,col_list(j)) dtip.y(k,col_list(j))],...
            'String',sprintf('X: %s\nY: %s',...
            num2str(dtip.x(k,col_list(j))),num2str(dtip.y(k,col_list(j)))));
    end
end

