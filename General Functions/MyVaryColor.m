function ColorSet = MyVaryColor(NumberOfSets, usrColormap)
    if nargin<2
        currentMap = jet;
%         currentMap = My_colorMap();
    else
        currentMap = usrColormap;
    end
    
    ColorMapSize = size(currentMap, 1);
    if NumberOfSets < ColorMapSize
        select_from_colormap = linspace(1,ColorMapSize,NumberOfSets);
        ColorSet(:,1) = interp1(1:ColorMapSize,currentMap(:,1),select_from_colormap);
        ColorSet(:,2) = interp1(1:ColorMapSize,currentMap(:,2),select_from_colormap);
        ColorSet(:,3) = interp1(1:ColorMapSize,currentMap(:,3),select_from_colormap);
    elseif NumberOfSets == ColorMapSize
        ColorSet = currentMap;
    else
        scale = linspace(1,NumberOfSets, ColorMapSize);
        ColorSet(:,1) = interp1(scale,currentMap(:,1),1:NumberOfSets);
        ColorSet(:,2) = interp1(scale,currentMap(:,2),1:NumberOfSets);
        ColorSet(:,3) = interp1(scale,currentMap(:,3),1:NumberOfSets);
    end

end