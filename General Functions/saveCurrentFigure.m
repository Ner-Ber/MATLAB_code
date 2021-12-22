function saveCurrentFigure(directory,Name)
    
    FigHandle = gcf;
    if nargin<2
        A = get(FigHandle,'children');
        AA = A(arrayfun(@(x) isa(x,'matlab.graphics.axis.Axes'), A));
        FigTitle = AA.Title.String;
        FigTitle = strrep(FigTitle,'{','');
        FigTitle = strrep(FigTitle,'}','');
        FigTitle = strrep(FigTitle,'^','');
        FigTitle = strrep(FigTitle,'\','');
        FigTitle = strrep(FigTitle,'[','');
        FigTitle = strrep(FigTitle,']','');
        FigTitle = strrep(FigTitle,',','');
        FigTitle = strrep(FigTitle,'.','');
        FigTitle = strrep(FigTitle,'*','x');
        Name = FigTitle;
    end
    set(FigHandle,'Color','w');
    savefig(FigHandle, fullfile(directory, [Name, '.fig']));
    export_fig(FigHandle,fullfile(directory, [Name, '.png']),'-png','-r300') 
end