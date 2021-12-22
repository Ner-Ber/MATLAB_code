close all
for i=find(relevantCells)
    %% make the new directory
    motherFolder = 'C:\Users\NeriB\Google Drive\JAY_lab\PROGRESS\eventsComparisons';
    Ev=i;
    ExpHour = PhediStructCellUpdated{i}.ExperimentData.ExpHour;
    ExpHour(regexp(ExpHour,'-'))=[];
    ExpDate = PhediStructCellUpdated{i}.ExperimentData.ExpDate;
    ExpDate(regexp(ExpDate,'-'))=[];
    fileName = [num2str(Ev),'_',num2str(ExpHour),'_',num2str(ExpDate)];
    
%     figure; IDT_PlotRowOverTime(PhediStructCellUpdated{i}.BigPicRotStruct);
%     newDir = fullfile(motherFolder,'BigPic');
%     mkdir(newDir);
%     saveas(gcf,fullfile(newDir,[fileName,'.fig']));
    
    figure; hold on;
    plot(PhediStructCellUpdated{i}.SgData.x_mins_x_tips(:,8),1e-3*(PhediStructCellUpdated{i}.SgData.Uxx(:,8)-mean(PhediStructCellUpdated{i}.SgData.Uxx(1:100,8))),'.-','Color',rgb('Pink'));
    plot(PhediStructCellUpdated{i}.solAtSG.x,PhediStructCellUpdated{i}.solAtSG.Uxx,'r');
    plot(PhediStructCell{i}.solAtSG.x,PhediStructCell{i}.solAtSG.Uxx,'b');
    title({PhediStructCell{i}.BigPicRotStruct.details,...
        ['\Gamma_{new}=',num2str(PhediStructCellUpdated{i}.solAtSG.Gamma),'  \Gamma_{old}=',num2str(PhediStructCell{i}.solAtSG.Gamma)]});
    xlim([-0.05 0.05])
    newDir = fullfile(motherFolder,'sg');
    mkdir(newDir);
    saveas(gcf,fullfile(newDir,[fileName,'.png']));
    
%     Movie_phedi_plotWithSg(PhediStructCellUpdated{i},'11',{'ROT'});
%     newDir = fullfile(motherFolder,'ROT');
%     mkdir(newDir);
%     saveas(gcf,fullfile(newDir,[fileName,'.fig']));
%     
%     Movie_phedi_plotWithSg(PhediStructCellUpdated{i},'11',{'locSpace'});
%     xlim([-0.05 0.05]);
%     newDir = fullfile(motherFolder,'locSpace');
%     mkdir(newDir);
%     saveas(gcf,fullfile(newDir,[fileName,'.fig']));
%     
%     Movie_phedi_plotWithSg(PhediStructCellUpdated{i},'11',{'velSpace'});
%     xlim([-0.05 0.05])
%     newDir = fullfile(motherFolder,'velSpace');
%     mkdir(newDir);
%     saveas(gcf,fullfile(newDir,[fileName,'.fig']));
    
    close all
end