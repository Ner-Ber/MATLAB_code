function hourStr = exp_getHourFromDataStructCell(DataStructCell,ExpNum)
    L = length(DataStructCell);
    l = 0;
    hourStr = [];
    while isempty(hourStr) && l<L
        try
            hourStr = strrep(DataStructCell{l+1}.ExperimentData.ExpHour,'-','');
            l = l+1;
        catch 
            l = l+1;
            continue
        end
    end
    if isempty(hourStr)
        hourStr = ['ExpNum',num2str(ExpNum)];
    end
    
end