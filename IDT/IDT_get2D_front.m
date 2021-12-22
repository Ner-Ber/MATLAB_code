function [uniqTimes, uniqLocs] = IDT_get2D_front(exp_dir,eventNum,varargin)

[pre_time,intervalTime,post_time,smt,numOfLines] =...
    setDefaults4function(varargin,3e-3,1e-6,3e-3,0,8);
%% iterate on lines
allFrontTimes = cell(1,numOfLines);
allFrontLocs = cell(1,numOfLines);
for lineNum = 1:numOfLines
    %--- get ims structure per row
    ims = IDT_imagesWithData(exp_dir,eventNum,pre_time,intervalTime,post_time,smt,lineNum);
    %--- find front per row
    frontStructure = IDT_FindCfFromRowOverTime(ims.lines);
%     frontIntepStruct = IDT_interpolateFrontBetweenSteps(frontStructure);
    frontIntepStruct = IDT_interpolateFrontOnFrames(frontStructure);
    allFrontTimes{lineNum} = frontIntepStruct.tq; %+ims.t(1)*ims.fps-1;    % time vector for front in units of frames
    allFrontLocs{lineNum} = frontIntepStruct.xq;
end

%% organize the data
ALL_times = cell2mat(allFrontTimes);
uniqTimes = unique(ALL_times);
uniqLocs = nan(numOfLines,length(uniqTimes));
for lineNum = 1:numOfLines
    Lia = ismember(uniqTimes,allFrontTimes{lineNum});
    uniqLocs(lineNum,Lia) = allFrontLocs{lineNum};
end


end