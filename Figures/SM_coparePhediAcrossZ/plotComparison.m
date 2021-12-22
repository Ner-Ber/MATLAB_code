%% compare phedi behavior across lines (z axis)
% purpose of this script is to create figures comparing phedi behavior
% across z axis.
% the data of each side of the block (large z and small z) is created using
% the follwoing script:
% batchGetDataStruct_specificLine.m

%% load regular avg data
% avgDatPath = 'F:\NeriFricsData\Frics\2018-10-10\allPhediStructures_17186.mat';
avgDatPath = 'F:\NeriFricsData\Frics\2018-10-10\allPhediStructures_17186_OLD.mat';
% avgDatPath = 'E:\Frics\2018-10-10\allPhediStructures_17186.mat';
load(avgDatPath);
relevantEvents = find(cellfun(@isRelevantEvent,PhediStructCell));
relevantEvents = relevantEvents(relevantEvents~=1);
relevantEvents = [4,7,16]; %--- events that are good examples for sync along z
for iEv = 1:length(relevantEvents)
    PhediStructCell{relevantEvents(iEv)}.PhediDataSG = PhediStructCell{relevantEvents(iEv)}.PhediData;
    PhediStructCell{relevantEvents(iEv)}.PhediData = PhediStructCell{relevantEvents(iEv)}.PhediDataSimpSmth;
end

Cf_vec = [];
for iEv = 1:length(relevantEvents)
    Cf_vec(iEv) = PhediStructCell{relevantEvents(iEv)}.PhediData.Cf;
end
relevantEvents(Cf_vec>1260) = [];
Cf_vec(Cf_vec>1260) = [];

%% load and plot the difference on z axis
z_diff_PhediStruct_path = 'F:\NeriFricsData\Frics\2018-10-10\17-18-6_sepLines';
z_diff_PhediStruct_nameTemplate = 'PhediStructures_multipLines_exp_17-18-6_ev_%d.mat';
for Ev = relevantEvents(:)'
% for Ev = 16
    z_diff_struct = load(fullfile(z_diff_PhediStruct_path, sprintf(z_diff_PhediStruct_nameTemplate,Ev)));   % load the specific file
    
    %-- if one side of the z axis is missing pass to the next event
    if isfield(z_diff_struct.PhediStructMultipLines.line_1_2_3,'theWorkspace') || isfield(z_diff_struct.PhediStructMultipLines.line_6_7_8,'theWorkspace')
        continue
    end
    
%     z_diff_struct.PhediStructMultipLines.line_1_2_3.PhediData = z_diff_struct.PhediStructMultipLines.line_1_2_3.PhediDataSimpSmth;
%     z_diff_struct.PhediStructMultipLines.line_6_7_8.PhediData = z_diff_struct.PhediStructMultipLines.line_6_7_8.PhediDataSimpSmth;
%     z_diff_struct.PhediStructMultipLines.line_1_2_3.PhediData = z_diff_struct.PhediStructMultipLines.line_1_2_3.PhediDataPreSmth;
%     z_diff_struct.PhediStructMultipLines.line_6_7_8.PhediData = z_diff_struct.PhediStructMultipLines.line_6_7_8.PhediDataPreSmth;
    
    [AvgPhediStruct_1_2_3, AvgPhediStruct_6_7_8, AvgPhediStruct_all] = plot_z_diff_comparison(z_diff_struct,PhediStructCell, Ev, Cf_vec(relevantEvents==Ev));

        
    Movie_phedi_plotWithSg(z_diff_struct.PhediStructMultipLines.line_1_2_3,'00',{'ROT','velSpace'})
    Movie_phedi_plotWithSg(z_diff_struct.PhediStructMultipLines.line_6_7_8,'00',{'ROT','velSpace'})
    Movie_phedi_plotWithSg(PhediStructCell{Ev},'00',{'ROT','velSpace'})
    
end