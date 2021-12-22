function DataStruct = Movie_createDataStructOfBigSmallPics(experNum, eventNum,Param)
    % [DataStruct] = Movie_createDataStructOfBigSmallPics(experNum, eventNum,varargin)
    %
    % Movie_createDataStructOfBigSmallPics will read the images from both
    % camers and create the relevant data structures.
    %
    % Insert the optional parameters as a structure with the fllowing fileds
    % (presented with recomended values):
    %
    %     'experNum',             1;...
    %     'eventNum',             3;...
    %     'lineNum',              'all';...
    %     'BlurROT',              1;...
    %     'UsrRes',               nan;...
    %     'numOOM',               3;...
    %     'slopePenalty',         0;...
    %     'expandLowBound',       -1;...
    %     'expandHighBound',      0;...
    %     'followMethod',         'RMS_minimization';...
    %     'phediRelative2first',  0;...
    %     'coarseShift',          0;...
    %     'preRowsTime',          3e-3;...
    %     'postRowsTime',         3e-3;...
    %     'prePhediTime',         1.5e-4;...
    %     'postPhediTime',        3e-4;...
    %     'smooth4Phedi',         0.001;...
    %     'PhediSafetyDist',      30;...
    %     'QuickAndDirty',        0;...
    %     'scratchWidth',         100*1e-6;...
    %     'scratchRegionMeters',  [0.07 0.09];...
    %     'defPhotoLoc',          0.08;...
    %     'Cr',                   1255;...
    %     'CohsvModl_h',          1e-7;...
    %     'skipPrecursor',        1;...
    %     'frontFallDetermn',     2
   
        %% general:
        %--- experiment path
        exper = my_dir;
        exp_dir = exper{experNum};
        
        %--- create the data struct
        DataStruct = struct;
        DataStruct.Param = Param;
        %% Camera meta data
        CamMetaStruct = CameraMetaAllCams(exp_dir);
        %--- update data structure
        DataStruct.CamMeta = CamMetaStruct;
        
        %% Get row over time
        [RowOverTime, timeCount] = phantom_getRowOverTime_and_timeCount(exp_dir, eventNum,Param.preRowsTime, Param.postRowsTime, Param.lineNum);
        %     RowOverTime = phantom_getRowOverTime(exp_dir, eventNum, CamMetaStruct, Param.preRowsTime, Param.postRowsTime, Param.lineNum);
        %--- blur rows of ROT if needed
        RowOverTimePreBlur = RowOverTime;
        if Param.BlurROT>1 && mod(Param.BlurROT,1)==0 && isnumeric(Param.BlurROT)
            g = createGaussianAproxCostume(1,Param.BlurROT);
            RowOverTime = conv2(RowOverTime,g,'same');
        end
        %% get resolution
        if isnan(Param.UsrRes)
            expMeasures = expMeasureProperties(exp_dir);
%             res = phantom_getResFromScrathces(RowOverTime);
            res = expMeasures.resSmallPic;
        else
            res = Param.UsrRes;
        end
        
        %% create Experiment Data structure
        [~,ExpDate,~] = fileparts(cd);
        ExpHour = exp_dir;
        
        ExperimentData = struct(...
            'experNum',experNum,...
            'eventNum',eventNum,...
            'lineNum',Param.lineNum,...
            'res',res,...
            'numOOM',Param.numOOM,...
            'preRowsTime',Param.preRowsTime,...
            'postRowsTime',Param.postRowsTime,...
            'prePhediTime',Param.prePhediTime,...
            'postPhediTime',Param.postPhediTime,...
            'scratchRegionMeters',Param.scratchRegionMeters,...
            'ExpDate',ExpDate,...
            'ExpHour',ExpHour...
            );
        %--- update the data struct
        exp_details=expDetailsRead(exp_dir);
        f = fieldnames(exp_details);
        for i = 1:length(f)
            ExperimentData.(f{i}) = exp_details.(f{i});
        end
        DataStruct.ExperimentData = ExperimentData;
        
        
        %% follow and analyze asperities
        disp('Definning and following asperities...');
        [vallyPairs,~,maxTrajectories,CM_Trajectories,skewnessMat,varMat,center_of_mass, MarkedMatCell]...
            = Movie_followAsperities(RowOverTime);
        
        %% create a normalized RowOverTime foe improved phedi-following
        RowOverTime_normalized = ROT_normalizeIntesByAsperity(RowOverTime,MarkedMatCell);
        
        %%
        fps = CamMetaStruct.FrameRate;        % frames/s
        postFrames = round(Param.postRowsTime*fps);
        preFrames = round(Param.preRowsTime*fps);
        %--- regulate number of post/pre frames
        if postFrames>CamMetaStruct.PostIms; postFrames=CamMetaStruct.PostIms; end
        if preFrames>(CamMetaStruct.NumIms-CamMetaStruct.PostIms); postFrames=CamMetaStruct.PostIms; end
        
        spatialVec = (0:(size(RowOverTime,2)-1))/res;
        
        %% save asperity data structure
        disp('saving asperity data...');
        analyzeAsperityStruct = struct(...
            'RowOverTime',RowOverTime,...
            'RowOverTimePreBlur', RowOverTimePreBlur,...
            'RowOverTime_normalized',RowOverTime_normalized,...
            'vallyPairs',vallyPairs,...
            'maxTrajectories',{maxTrajectories},...
            'CM_Trajectories',{CM_Trajectories},...
            'skewness',skewnessMat,...
            'var',varMat,...
            'center_of_mass',center_of_mass,...
            'spatialVec',spatialVec,...
            'CameraResolution',res,...
            'timeCount',timeCount,...
            'preFrames',preFrames,...
            'postFrames',postFrames...
            );
        
        %--- update the data struct
        DataStruct.AsperityData = analyzeAsperityStruct;
        
        %% get Big Picture data
        SubFolders = my_dir(exp_dir);
        if (sum(strcmpi(SubFolders,'Ph'))+sum(strcmpi(SubFolders,'PhBig')))==2
            BigPicRowOverTimeStruct = phantom_BigPicStruct(exp_dir, eventNum,Param.preRowsTime,1e-10,Param.postRowsTime,1,1,'all',Param.frontFallDetermn);
        else
            BigPicRowOverTimeStruct = IDT_PlotRowOverTime(exp_dir,eventNum,Param.preRowsTime,0,Param.postRowsTime,0,1:CamMetaStruct.ImageHeight,[],[],0);
        end
        %--- update the data struct
        DataStruct.BigPicRotStruct = BigPicRowOverTimeStruct;
        
    end