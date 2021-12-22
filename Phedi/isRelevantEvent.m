function Logic = isRelevantEvent(DataStruct)
    if isfield(DataStruct,'theWorkspace') Logic=false; else Logic=~ischar(DataStruct.PhediData); end
end