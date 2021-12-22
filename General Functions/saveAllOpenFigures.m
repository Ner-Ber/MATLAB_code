function saveAllOpenFigures(directory, varargin)
    
    [useFigName] = setDefaults4function(varargin,0);
    
    % directory = 'G:\Frics\2018-10-10\Exp7Figs';   % Your destination folder
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
        FigHandle = FigList(iFig);
        if useFigName
            FigTitle = FigHandle.Name;
        else
            A = get(FigHandle,'children');
            AA = A(arrayfun(@(x) isa(x,'matlab.graphics.axis.Axes'), A));
            FigTitle = AA.Title.String;
            FigTitle = FigTitle';
            FigTitle = FigTitle(:)';
    %         FigTitle = FigTitle(2,:);
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
        export_fig(FigHandle,fullfile(directory, [FigTitle, '.png']),'-png','-r150')
%         savefig(FigHandle, fullfile(directory, [num2str(FigList(iFig).Number), '.fig']));
%         export_fig(FigHandle,fullfile(directory, [num2str(FigList(iFig).Number), '.png']),'-png','-r300')
        %         T = [FigTitle{1},' ',FigTitle{2}];
        %         savefig(FigHandle, fullfile(directory, [T, '.fig']));
        %         export_fig(FigHandle,fullfile(directory, [T, '.png']),'-png','-r300')
    end
    
end