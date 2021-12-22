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
pderect([-1.4954909819639282 1.4984969939879762 -0.85220440881763526 -0.99649298597194402],'R1');
set(findobj(get(pde_fig,'Children'),'Tag','PDEEval'),'String','R1')

% Boundary conditions:
pdetool('changemode',0)
pdesetbd(4,...
'dir',...
2,...
str2mat('1','0','0','1'),...
str2mat('0','0'))
pdesetbd(3,...
'neu',...
2,...
str2mat('0','0','0','0'),...
str2mat('0','0'))
pdesetbd(2,...
'dir',...
2,...
str2mat('1','0','0','1'),...
str2mat('-0.003','0'))
pdesetbd(1,...
'neu',...
2,...
str2mat('0','0','0','0'),...
str2mat('0','0'))

% Mesh generation:
setappdata(pde_fig,'Hgrad',1.3);
setappdata(pde_fig,'refinemethod','regular');
setappdata(pde_fig,'jiggle',char('on','mean',''));
pdetool('initmesh')
pdetool('refine')
pdetool('refine')
pdetool('refine')
pdetool('refine')

% PDE coefficients:
pdeseteq(1,...
str2mat('2*((3E8)./(2*(1+(0.3))))+(2*((3E8)./(2*(1+(0.3)))).*(0.3)./(1-(0.3)))','0','(3E8)./(2*(1+(0.3)))','0','(3E8)./(2*(1+(0.3)))','2*((3E8)./(2*(1+(0.3)))).*(0.3)./(1-(0.3))','0','(3E8)./(2*(1+(0.3)))','0','2*((3E8)./(2*(1+(0.3))))+(2*((3E8)./(2*(1+(0.3)))).*(0.3)./(1-(0.3)))'),...
str2mat('0.0','0.0','0.0','0.0'),...
str2mat('0.0','0.0'),...
str2mat('1.0','0','0','1.0'),...
'0:10',...
'0.0',...
'0.0',...
'[0 100]')
setappdata(pde_fig,'currparam',...
['3E8';...
'0.3';...
'0.0';...
'0.0';...
'1.0'])

% Solve parameters:
setappdata(pde_fig,'solveparam',...
str2mat('0','11520','10','pdeadworst',...
'0.5','longest','0','1E-4','','fixed','Inf'))

% Plotflags and user data strings:
setappdata(pde_fig,'plotflags',[8 1 8 1 1 1 6 1 0 0 0 1 1 0 0 0 0 1]);
setappdata(pde_fig,'colstring','');
setappdata(pde_fig,'arrowstring','');
setappdata(pde_fig,'deformstring','');
setappdata(pde_fig,'heightstring','');

% Solve PDE:
pdetool('solve')