function operate_on_all_figs(command, varargin)
% change some property on all or some of the fgures. figs2change (2nd
% argument) can be either 'all' or a vector containing figure numbers.
% command should be in string format

[figs2change] = setDefaults4function(varargin,'all');

all_figs = findobj(0, 'type', 'figure');
if ischar(figs2change)
    fig_list = all_figs;
else
    fig_list = figs2change;
end

InterFig = intersect(all_figs, fig_list);
for f = 1:length(InterFig)
    set(0, 'CurrentFigure', InterFig(f));
    eval(command);
end


end