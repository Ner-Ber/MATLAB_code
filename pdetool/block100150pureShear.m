function pdemodel
[pde_fig,ax]=pdeinit;
pdetool('appl_cb',3);
set(ax,'DataAspectRatio',[1 1 1]);
set(ax,'PlotBoxAspectRatio',[1.5 1 1]);
set(ax,'XLim',[-1.5 1.5]);
set(ax,'YLim',[-1 1]);
set(ax,'XTickMode','auto');
set(ax,'YTickMode','auto');

% Geometry description:
pderect([-0.99552117263843654 0.49959283387622122 -0.8135179153094465 0.19543973941368065],'R1');
set(findobj(get(pde_fig,'Children'),'Tag','PDEEval'),'String','R1')

% Boundary conditions:
pdetool('changemode',0)
pdesetbd(4,...
'neu',...
2,...
str2mat('0','0','0','0'),...
str2mat('0','0'))
pdesetbd(3,...
'dir',...
2,...
str2mat('1','0','0','1'),...
str2mat('0','0'))
pdesetbd(2,...
'neu',...
2,...
str2mat('0','0','0','0'),...
str2mat('0','0'))
pdesetbd(1,...
'dir',...
2,...
str2mat('1','0','0','1'),...
str2mat('1','0'))

% Mesh generation:
setappdata(pde_fig,'Hgrad',1.3);
setappdata(pde_fig,'refinemethod','regular');
setappdata(pde_fig,'jiggle',char('on','mean',''));
pdetool('initmesh')
pdetool('refine')

% PDE coefficients:
pdeseteq(1,...
str2mat('2*((3E8)./(2*(1+(0.33))))+(2*((3E8)./(2*(1+(0.33)))).*(0.33)./(1-(0.33)))','0','(3E8)./(2*(1+(0.33)))','0','(3E8)./(2*(1+(0.33)))','2*((3E8)./(2*(1+(0.33)))).*(0.33)./(1-(0.33))','0','(3E8)./(2*(1+(0.33)))','0','2*((3E8)./(2*(1+(0.33))))+(2*((3E8)./(2*(1+(0.33)))).*(0.33)./(1-(0.33)))'),...
str2mat('0.0','0.0','0.0','0.0'),...
str2mat('0.0','0.0'),...
str2mat('1.0','0','0','1.0'),...
'0:10',...
'0.0',...
'0.0',...
'[0 100]')
setappdata(pde_fig,'currparam',...
['3E8 ';...
'0.33';...
'0.0 ';...
'0.0 ';...
'1.0 '])

% Solve parameters:
setappdata(pde_fig,'solveparam',...
str2mat('0','1392','10','pdeadworst',...
'0.5','longest','0','1E-4','','fixed','Inf'))

% Plotflags and user data strings:
setappdata(pde_fig,'plotflags',[12 1 1 1 1 1 7 1 0 0 0 1 1 0 0 1 0 1]);
setappdata(pde_fig,'colstring','');
setappdata(pde_fig,'arrowstring','');
setappdata(pde_fig,'deformstring','');
setappdata(pde_fig,'heightstring','');

% Solve PDE:
pdetool('solve')
