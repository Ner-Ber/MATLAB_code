%% set request for data
% lineNum_cell = {1,2,3,4,5,6,7,8};
% lineNum_cell = {1:3,6:8, 'all'};
lineNum_cell = {1:4,5:8, 'all'};
date = '2018-10-10';
cd(fullfile('F:\NeriFricsData\Frics\',date));

%% experiments to read
exp_name = {'17-18-6'};
exp_num = nan(1,length(exp_name));
for i=1:length(exp_name)
    exp_num(i) = find(strcmp(exp_name{i},my_dir));
end
% events = {[4     6     7     8     9    12    14    15    16    17    20    21    22    24    25    30]};
% events = {[12    16    17]};
events = {[7,4,12,16,17,22,24]};
signalSmooth4PhediTrack = 1;
blurROT = 1;
%% get new data
for i = 1:length(exp_num)
    ExperNum = exp_num(i);
    for EvNum = events{i}
        PhediStructMultipLines = analyze_batchAnalyzePhedis_sepLines(ExperNum, EvNum, lineNum_cell, signalSmooth4PhediTrack, blurROT);
        
        saveFolder = 'F:\NeriFricsData\Frics\2018-10-10\17-18-6_sepLines';
        fileName = sprintf('PhediStructures_multipLines_2_exp_%s_ev_%d.mat',exp_name{i},EvNum);
        name = fullfile(saveFolder, fileName);
        save(name,'PhediStructMultipLines','-v7.3');
    end
end