close all
n=4;
m=5;
%% plot sg
SGfig = figure;
ha = my_tight_subplot(n,m,[.07 .02],[.01 .1],[.01 .01],0);
for idx=1:nnz(relevantCells)
    F = find(relevantCells);
    i = F(idx);
    
    axes(ha(idx));
    hold on;
    plot(PhediStructCellUpdated{i}.SgData.x_mins_x_tips(:,8),1e-3*(PhediStructCellUpdated{i}.SgData.Uxx(:,8)-mean(PhediStructCellUpdated{i}.SgData.Uxx(1:100,8))),'.-','Color',rgb('Pink'));
    plot(PhediStructCellUpdated{i}.solAtSG.x,PhediStructCellUpdated{i}.solAtSG.Uxx,'r');
    title(['Ev=',num2str(i),'  \Gamma_{new}=',num2str(PhediStructCellUpdated{i}.solAtSG.Gamma),'  Cf=',num2str(PhediStructCellUpdated{i}.PhediData.Cf)],'FontSize',9);
    xlim([-0.05 0.05])
    hold off
end

annotation('textbox', [0 0.9 1 0.1], ...
    'String', PhediStructCellUpdated{i}.BigPicRotStruct.details(1:23), ...
    'FontWeight','bold',...
    'FontSize',10,...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')

%--- save
set(SGfig, 'Color', 'w');
direct = 'C:\Users\NeriB\Google Drive\JAY_lab\PROGRESS\eventsComparisons\sg';
name = PhediStructCellUpdated{i}.BigPicRotStruct.details(1:23);
fullName = fullfile(direct,name);
export_fig(SGfig,fullName,'-png','-r600')
savefig(SGfig,[fullName,'.fig']);


%% plot ROT
ROTfig = figure;
ha = my_tight_subplot(n,m,[.07 .02],[.01 .1],[.01 .01],0);
for idx=1:nnz(relevantCells)
    F = find(relevantCells);
    i = F(idx);
    
%     axes(ha(idx));
    hold on;
    Movie_phedi_plotWithSg(PhediStructCellUpdated{i},'11',{'ROT'},[],[],ha(idx));
    ha(idx).XLabel.String = '';
    ha(idx).YLabel.String = '';
    title(['Ev=',num2str(i),'  Cf=',num2str(PhediStructCellUpdated{i}.PhediData.Cf)],'FontSize',9);
    hold off
end

annotation('textbox', [0 0.9 1 0.1], ...
    'String', PhediStructCellUpdated{i}.BigPicRotStruct.details(1:23), ...
    'FontWeight','bold',...
    'FontSize',10,...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')

%--- save
set(ROTfig, 'Color', 'w');
direct = 'C:\Users\NeriB\Google Drive\JAY_lab\PROGRESS\eventsComparisons\ROT';
name = PhediStructCellUpdated{i}.BigPicRotStruct.details(1:23);
fullName = fullfile(direct,name);
export_fig(ROTfig,fullName,'-png','-r600')
savefig(ROTfig,[fullName,'.fig']);


%% plot loc
locFig = figure;
ha = my_tight_subplot(n,m,[.07 .02],[.01 .1],[.01 .01],0);
for idx=1:nnz(relevantCells)
    F = find(relevantCells);
    i = F(idx);
    
%     axes(ha(idx));
    hold on;
    Movie_phedi_plotWithSg(PhediStructCellUpdated{i},'11',{'locSpace'},[],[],ha(idx));
    ha(idx).XLabel.String = '';
    ha(idx).YLabel.String = '';
    title(['Ev=',num2str(i),'  Cf=',num2str(PhediStructCellUpdated{i}.PhediData.Cf)],'FontSize',9);
    xlim([-0.05 0.05]);
    hold off
end

annotation('textbox', [0 0.9 1 0.1], ...
    'String', PhediStructCellUpdated{i}.BigPicRotStruct.details(1:23), ...
    'FontWeight','bold',...
    'FontSize',10,...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')

%--- save
set(locFig, 'Color', 'w');
direct = 'C:\Users\NeriB\Google Drive\JAY_lab\PROGRESS\eventsComparisons\locSpace';
name = PhediStructCellUpdated{i}.BigPicRotStruct.details(1:23);
fullName = fullfile(direct,name);
export_fig(locFig,fullName,'-png','-r600')
savefig(locFig,[fullName,'.fig']);


%% plot vel
velFig = figure;
ha = my_tight_subplot(n,m,[.07 .02],[.01 .1],[.01 .01],0);
for idx=1:nnz(relevantCells)
    F = find(relevantCells);
    i = F(idx);
    
%     axes(ha(idx));
    hold on;
    Movie_phedi_plotWithSg(PhediStructCellUpdated{i},'11',{'velSpace'},[],[],ha(idx));
    ha(idx).XLabel.String = '';
    ha(idx).YLabel.String = '';
    title(['Ev=',num2str(i),'  Cf=',num2str(PhediStructCellUpdated{i}.PhediData.Cf)],'FontSize',9);
    xlim([-0.05 0.05])
    hold off
end

annotation('textbox', [0 0.9 1 0.1], ...
    'String', PhediStructCellUpdated{i}.BigPicRotStruct.details(1:23), ...
    'FontWeight','bold',...
    'FontSize',10,...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')

%--- save
set(velFig, 'Color', 'w');
direct = 'C:\Users\NeriB\Google Drive\JAY_lab\PROGRESS\eventsComparisons\velSpace';
name = PhediStructCellUpdated{i}.BigPicRotStruct.details(1:23);
fullName = fullfile(direct,name);
export_fig(velFig,fullName,'-png','-r600')
savefig(velFig,[fullName,'.fig']);

%% plot BigPicROT
BigPicFig = figure;
ha = my_tight_subplot(n,m,[.07 .02],[.01 .1],[.01 .01],0);
for idx=1:nnz(relevantCells)
    F = find(relevantCells);
    i = F(idx);
    
%     axes(ha(idx));
    hold on;
    IDT_PlotRowOverTime(PhediStructCellUpdated{i}.BigPicRotStruct,[],[],[],[],[],[],[],ha(idx));
    colorbar off;
    ha(idx).XLabel.String = '';
    ha(idx).YLabel.String = '';
    title(['Ev=',num2str(i),'  Cf=',num2str(PhediStructCellUpdated{i}.PhediData.Cf)],'FontSize',9);
    ylim(([-5e-4 5e-4]+mean(PhediStructCellUpdated{i}.PhediData.t_tips)))
    xlim([0.06 0.1]);
    hold off
end

annotation('textbox', [0 0.9 1 0.1], ...
    'String', PhediStructCellUpdated{i}.BigPicRotStruct.details(1:23), ...
    'FontWeight','bold',...
    'FontSize',10,...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')

%--- save
set(BigPicFig, 'Color', 'w');
direct = 'C:\Users\NeriB\Google Drive\JAY_lab\PROGRESS\eventsComparisons\BigPic';
name = PhediStructCellUpdated{i}.BigPicRotStruct.details(1:23);
fullName = fullfile(direct,name);
export_fig(BigPicFig,fullName,'-png','-r600')
savefig(BigPicFig,[fullName,'.fig']);