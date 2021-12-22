%% complete first event
ExpNum = 7;
PhediStructCell = analyze_batchAnalyzePhedis(ExpNum,10);
%-- save
hourStr = exp_getHourFromDataStructCell(PhediStructCell,ExpNum);
name=['allPhediStructures_',hourStr,'.mat'];
save(name,'PhediStructCell','-v7.3');

%% get new data
for i= [4:6,8]
    PhediStructCell = analyze_batchAnalyzePhedis(i);
    hourStr = exp_getHourFromDataStructCell(PhediStructCell,i);
    name=['allPhediStructures_',hourStr,'.mat'];
    save(name,'PhediStructCell','-v7.3');
end