%% set request for data
lineNum_vec = {1,2,3,4,5,6,7,8};
date = '2018-10-10';
cd(fullfile('F:\NeriFricsData\Frics\',date));

%% experiments to read
exp_date = {'17-18-6'};
% exp_num = nan(1,length(exp_name));
% for i=1:length(exp_name)
%     exp_num(i) = find(strcmp(exp_name{i},my_dir));
% end
events = {[4     6     7     8     9    12    14    15    16    17    20    21    22    24    25    30]};

%% get new data


for i_dates = 1:length(exp_date)
    allPhediStructures_existing = load(strcat('allPhediStructures_',strrep(exp_date{i_dates},'-',''),'.mat'));
    allPhediStructures_existing = allPhediStructures_existing.PhediStructCell;
%     BigPicMultipLines_structs = cell(size(allPhediStructures_existing));
    for i_event = events{i_dates}
        %--- get original Params structure
        Param = allPhediStructures_existing{i_event}.Param;
        BigPicMultipLines = struct;
        for l_idx = 1:length(lineNum_vec)
            BigPicRowOverTimeStruct = phantom_BigPicStruct(exp_date{i_dates}, i_event,Param.preRowsTime,1e-10,Param.postRowsTime,1,1,...
                lineNum_vec{l_idx},Param.frontFallDetermn);      
            BigPicMultipLines.(strcat('line_',num2str(lineNum_vec{l_idx}))) = BigPicRowOverTimeStruct;
        end
        name = strcat('17-18-6_BigPics_byLines\BigPicMultipLines_Ev',num2str(i_event),'.mat');
        save(name,'BigPicMultipLines','-v7.3');
%         BigPicMultipLines_structs{i_event} = BigPicMultipLines;
    end
%     name = strcat('BigPicMultipLines_',strrep(exp_date{i_dates},'-',''),'.mat');
%     save(name,'BigPicMultipLines_structs','-v7.3');
end
