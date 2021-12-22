function close_all_except(figs2keep)
% close_all_except will close all figures except those specified by their
% number in 'figs2keep'

all_figs = findobj(0, 'type', 'figure');
delete(setdiff(all_figs, figs2keep));

end