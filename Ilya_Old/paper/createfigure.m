function createfigure(xdata1, ydata1, zdata1, X1, Y1, Y2)
%CREATEFIGURE(XDATA1,YDATA1,ZDATA1,CDATA1,X1,Y1,Y2)
%  XDATA1:  surface xdata
%  YDATA1:  surface ydata
%  ZDATA1:  surface zdata
%  CDATA1:  surface cdata
%  X1:  vector of x data
%  Y1:  vector of y data
%  Y2:  vector of y data

%  Auto-generated by MATLAB on 10-Jun-2013 12:49:24

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.0806451612903226 0.435714285714286 0.372503840245776 0.485714285714285],...
    'CLim',[0.6 1.1]);
%% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0 200.158478605388]);
%% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[-0.299325428903103 -0.098015397787094]);
%% Uncomment the following line to preserve the Z-limits of the axes
% zlim(axes1,[0.671264489161687 1.20614035087719]);
grid(axes1,'on');
hold(axes1,'all');

% Create mesh
mesh(xdata1,ydata1,zdata1,'Parent',axes1);
caxis([0.7 1.1]);
xlim([0 200]);

% Create xlabel
xlabel('X (mm)');

% Create title
title('2013-03-17 18-11-04  event - 9  Trigger=978.4224');

% Create colorbar
colorbar('peer',axes1);
% Resize the axes in order to prevent it from shrinking.
set(axes1,...
    'Position',[0.0806451612903226 0.435714285714286 0.372503840245776 0.485714285714285]);

% Create axes
axes2 = axes('Parent',figure1,...
    'Position',[0.559907834101382 0.528011204481793 0.335358898606797 0.400751063784292],...
    'FontSize',25);
box(axes2,'on');
hold(axes2,'all');

% Create plot
plot(X1,Y1,'Parent',axes2,'MarkerSize',20,'Marker','.');

% Create axes
axes3 = axes('Parent',figure1,...
    'Position',[0.559907834101382 0.0 0.335358898606797 0.400751063784292]);
box(axes3,'on');
hold(axes3,'all');

% Create plot
plot(X1,Y2,'Parent',axes3,'Marker','.');

% Create arrow
annotation(figure1,'arrow',[0.205069124423963 0.205837173579109],...
    [0.430616341030195 0.793960923623446],'HeadWidth',25,'HeadStyle','none',...
    'LineStyle','--',...
    'LineWidth',3);

