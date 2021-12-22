% saveFigures4Article(directory,figNums4Saving)
%
% saveFigures4Article will save figure in the amed directory in .emf and
% .fig formats, while applying a white background.
% If not specifically mentioned, saveFigures4Article will save all open
% figures. 
%

function saveFigures4Article(directory,figNums4Saving)
    
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    if nargin<2
        HandleIdxs = 1:length(FigList);
    else
        [~,HandleIdxs] = ismember(figNums4Saving,arrayfun(@(A) A.Number, FigList));
    end
    
    for iFig = HandleIdxs(:)'
        FigHandle = FigList(iFig);
        FigHandle.Color = 'w';
        FigTitle = FigHandle.Name;
        if isempty(FigTitle)
            FigTitle = ['subFig',num2str(FigHandle.Number)];
        end
        FigTitle = strrep(FigTitle,'{','');
        FigTitle = strrep(FigTitle,'}','');
        FigTitle = strrep(FigTitle,'^','');
        FigTitle = strrep(FigTitle,'\','');
        FigTitle = strrep(FigTitle,'[','');
        FigTitle = strrep(FigTitle,']','');
        FigTitle = strrep(FigTitle,'<','');
        FigTitle = strrep(FigTitle,'>','');
        FigTitle = strrep(FigTitle,'*','x');
        savefig(FigHandle, fullfile(directory, [FigTitle, '.fig']));
        saveas(FigHandle,fullfile(directory, [FigTitle, '.emf']));
        saveas(FigHandle,fullfile(directory, [FigTitle, '.svg']));
        saveas(FigHandle,fullfile(directory, [FigTitle, '.png']));
%         export_fig(FigHandle,fullfile(directory, [FigTitle, '.png']),'-png','-r300')
        %         T = [FigTitle{1},' ',FigTitle{2}];
        %         savefig(FigHandle, fullfile(directory, [T, '.fig']));
        %         export_fig(FigHandle,fullfile(directory, [T, '.png']),'-png','-r300')
    end
    
end