function my_createfigure(figNum)



% Create figure
%figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
set(0,'DefaultFigureWindowStyle','normal');
figure1=figure;
set(0,'DefaultFigureWindowStyle','docked');
set(gcf, 'units', 'centimeters', 'pos', [5 5 8 10])

Xsize=0.78;
Ysize=0.27;
X0=0.17;
Y0=0.66;

% ------------Create Upper axes
fig=my_get_axis(figNum(1));
%'YTick',[-0.4 0 0.4]
axes_tmp = axes('Parent',figure1,'XTick',fig.XTick,...
    'XTickLabel',{'','','','','','','','',''}, 'Position',[X0 Y0 Xsize Ysize]);
plot_Lines(axes_tmp, fig);
title(fig.title);

%------------ Create Middle axes
fig=my_get_axis(figNum(2));

%'YTick',[-1 -0.5 0]
axes_tmp = axes('Parent',figure1,'XTick',fig.XTick,...
    'XTickLabel',{'','','','','','','','',''}, 'Position',[X0 Y0-Ysize-0.01 Xsize Ysize]);
plot_Lines(axes_tmp, fig);
% ---------Create Bottom axes
fig=my_get_axis(figNum(3));
%'YTick',[-0.2 0 0.2]
axes_tmp = axes('Parent',figure1,'XTick',fig.XTick,...
    'XTickLabel',fig.XTick, 'Position',[X0 Y0-2*(Ysize+0.01) Xsize Ysize ]);
plot_Lines(axes_tmp, fig);

function plot_Lines(axes_tmp,fig)
xlim(axes_tmp,fig.XLim);
ylim(axes_tmp,fig.YLim);
box(axes_tmp,'on');
hold(axes_tmp,'all');
hold all;
ylabel(fig.yLabel);
xlabel(fig.xLabel);
%--if mat convert to cell
if( ~iscell(fig.x))
    
    if(length(fig.x(1,:))~=length(fig.y(1,:)))
        x=repmat(fig.x,1,length(fig.y(1,:)));
        y=fig.y;
        fig=rmfield(fig,'x');
        fig=rmfield(fig,'y');
    else
        x=fig.x;
        y=fig.y;
        fig=rmfield(fig,'x');
        fig=rmfield(fig,'y');
    end
    
    for j=1:length(x(1,:))
        fig.x{j}=(x(:,j));
        fig.y{j}=(y(:,j));
    end
    
end


for j=1:length(fig.x)
    a= plot(axes_tmp,fig.x{j},fig.y{j},'.-');
    set(a,'color',fig.color{j});
    set(a,'DisplayName',fig.DisplayName{j})
end
