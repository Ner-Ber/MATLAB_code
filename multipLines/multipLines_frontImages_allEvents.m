function multipLines_frontImages_allEvents(DAandREsidualsStruct, expHour)
    
    
    if ischar(DAandREsidualsStruct)
        DAandREsidualsStruct = load(DAandREsidualsStruct);
        DAandREsidualsStruct = DAandREsidualsStruct.myStruct;
    end
    
    thisExp = find(cellfun(@(A) strcmpi(A, ['2018-10-10 ',expHour]), extractfield(DAandREsidualsStruct, 'Name')));
    Cf_vec = DAandREsidualsStruct(thisExp).Cf_vec;
    Events = DAandREsidualsStruct(thisExp).Events;
    
    multipLinesFolder = strcat(expHour,'_BigPics_byLines');
    DIR = dir(multipLinesFolder);
    for d = 3:length(DIR)
        name = DIR(d).name;
        Ev = str2double(name(strfind(name,'_')+3:strfind(name,'.mat')-1));
        [multiLineData, SurfPlot, slicesFig] = multipLines_createImagesOfFront2D(fullfile(multipLinesFolder, name));
        figure(slicesFig);
        propString = sprintf('Ev=%d  C_f@phedis=%d',Ev, round(Cf_vec(Events==Ev)));
        suptitle(propString);
        SurfPlot.Parent.Title.String = propString;
        SurfPlot.Parent.Parent.Name = ['frontSurface_',propString([1:2,4:6])];
        slicesFig.Name = ['snapshots_',propString([1:2,4:6])];
    end
    
    
end