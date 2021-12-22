%% ***plot cohesive zone on Area intensity***
RowOverTimeStruct = PhediStruct.BigPicRotStruct;

refDist = [1,2.5,5,7.5,10,12]*(-1e-3);
cohesiveStructCell = {};
for i=1:length(refDist)
    cohesiveStructCell{i} = IDT_calcCohesiveFromROT(RowOverTimeStruct,[0.07 0.09],refDist(i));
%     cohesiveStructCell{i} = IDT_calcCohesiveFromROT_time2space(RowOverTimeStruct,[0.07 0.09],refDist(i));
    
    Colors = MyVaryColor(size(cohesiveStructCell{i}.renormalizedRows,1));
    %--- plot cohesive zone
    figure;
    subplot(1,2,1); hold on;
    L1 = [];
    LGND1 = {};
    for rowIdx = 1:size(cohesiveStructCell{i}.renormalizedRows,1)
        L1(rowIdx) = plot(cohesiveStructCell{i}.xVectorsMeters(rowIdx,:),...
            cohesiveStructCell{i}.renormalizedRows(rowIdx,:),...
            'LineWidth',1.2,'Color',Colors(rowIdx,:));
        plot([0 -cohesiveStructCell{i}.Xc(rowIdx)],[1 1]*(-0.63),'--',...
            'LineWidth',1.2,'Color',Colors(rowIdx,:));
        LGND1{rowIdx} = ['x_{tip}=',num2str(round(cohesiveStructCell{i}.x_tipsMeters(rowIdx),3)),...
            '  X_c=', num2str(round(cohesiveStructCell{i}.Xc(rowIdx),3))];
    end
    xlim([-0.02 0.01]); ylim([-2 0.2]);
    title(['cohesive with refDist=',num2str(refDist(i))]);
    legend(L1,LGND1,'Location','southeast');
    
    %--- plot dA zone
    subplot(1,2,2); hold on;
    L2 = [];
    LGND2 = {};
    for rowIdx = 1:size(cohesiveStructCell{i}.relevantFrames,1)
        L2(rowIdx) = plot(cohesiveStructCell{i}.xVectorsMeters(rowIdx,:),...
            cohesiveStructCell{i}.relevantFrames(rowIdx,:),...
            'LineWidth',1.2,'Color',Colors(rowIdx,:));
        plot([1 1]*refDist(i),[1 1-cohesiveStructCell{i}.dA(rowIdx)],'--',...
            'LineWidth',1.2,'Color',Colors(rowIdx,:));
        LGND2{rowIdx} = ['x_{tip}=',num2str(round(cohesiveStructCell{i}.x_tipsMeters(rowIdx),3)),...
            '  \deltaA=', num2str(round(cohesiveStructCell{i}.dA(rowIdx),3))];
    end
    xlim([-0.02 0.01]); ylim([0.7 1.1]);
    title(['dA with refDist=',num2str(refDist(i))]);
    legend(L2,LGND2,'Location','southeast');
end
