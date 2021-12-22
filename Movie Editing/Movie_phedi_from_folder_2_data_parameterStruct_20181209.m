function [DataStruct_WPhedis] = Movie_phedi_from_folder_2_data_parameterStruct_20181209(DataStruct)
    % [DataStruct] = Movie_phedi_from_folder_2_data_parameterStruct(experNum, eventNum,varargin)
    %
    % Movie_phedi_from_folder_2_data will read the images from the movie, then
    % claculate phedis and follow them end output cell arrays containig the
    % relevant data.
    %
    % Insert the optional parameters as a structure with the fllowing fileds
    % (presented with recomended values):
    %
    % 'experNum',             1
    % 'eventNum',             3
    % 'lineNum',              'all'
    % 'BlurROT',              20
    % 'UsrRes',               by phantom_getResFromScrathces function
    % 'numOOM',               3
    % 'phediRelative2first',  1
    % 'preRowsTime',          3e-3
    % 'postRowsTime',         3e-3
    % 'prePhediTime',         1.5e-4
    % 'postPhediTime',        3e-4
    % 'smooth4Phedi',         0.001
    % 'PhediSafetyDist',      30    (in Pix)
    % 'QuickAndDirty',        0
    % 'scratchWidth',         100*1e-6
    % 'scratchRegionMeters',  [0.07 0.09]
    % 'defPhotoLoc',          0.08
    % 'Cr',                   1255
    % 'CohsvModl_h',          1e-7
    % 'skipPrecursor',        1
    %
    %
    %
    %    03-12-2018 update:
    %       *) The function is now intended to work with a normalized
    %       RowOverTime produced by the function
    %       'ROT_normalizeIntesByAsperity'. It is produced in the previous
    %       function 'Movie_createDataStructOfBigSmallPics' and is inserted
    %       in the 'DataStruct.AsperityData' structure.
    %       *) Added option to cancel the coarse shift of row over time
    %
    
    
    
    try
        phediAnalysisTimer = tic;
        disp('analyzing phedis...');
        
        %% general and set variables
        DataStruct_WPhedis = DataStruct;
        Param = DataStruct.Param;
        CamMetaStruct = DataStruct.CamMeta;
        BigPicRotStruct = DataStruct.BigPicRotStruct;
        
        RowOverTime = DataStruct.AsperityData.RowOverTime;
        RowOverTime_normalized = DataStruct.AsperityData.RowOverTime_normalized;
        center_of_mass = DataStruct.AsperityData.center_of_mass;
        timeCount = DataStruct.AsperityData.timeCount;
        fps = CamMetaStruct.FrameRate;        % frames/s
        postFrames = round(Param.postRowsTime*fps);
        preFrames = round(Param.preRowsTime*fps);
        res = DataStruct.ExperimentData.res;
        
        %%
        %     %% general:
        %     %--- experiment path
        %     exper = my_dir;
        %     exp_dir = exper{experNum};
        %
        %     %--- create the data struct
        %     DataStruct = struct;
        %     DataStruct.Param = Param;
        %     %% Camera meta data
        %     CamMetaStruct = CameraMetaAllCams(exp_dir);
        %     %--- update data structure
        %     DataStruct.CamMeta = CamMetaStruct;
        %
        %     %% Get row over time
        %     [RowOverTime, timeCount] = phantom_getRowOverTime_and_timeCount(exp_dir, eventNum,Param.preRowsTime, Param.postRowsTime, Param.lineNum);
        %     %     RowOverTime = phantom_getRowOverTime(exp_dir, eventNum, CamMetaStruct, Param.preRowsTime, Param.postRowsTime, Param.lineNum);
        %     %--- blur rows of ROT if needed
        %     if Param.BlurROT>1 && mod(Param.BlurROT,1)==0 && isnumeric(Param.BlurROT)
        %         g = createGaussianAproxCostume(1,Param.BlurROT);
        %         RowOverTime = conv2(RowOverTime,g,'same');
        %     end
        %     %% get resolution
        %     if isnan(Param.UsrRes)
        %         res = phantom_getResFromScrathces(RowOverTime);
        %     else
        %         res = Param.UsrRes;
        %     end
        %
        %     %% create Experiment Data structure
        %     [~,ExpDate,~] = fileparts(cd);
        %     ExpHour = exp_dir;
        %
        %     ExperimentData = struct(...
        %         'experNum',experNum,...
        %         'eventNum',eventNum,...
        %         'lineNum',Param.lineNum,...
        %         'res',res,...
        %         'numOOM',Param.numOOM,...
        %         'preRowsTime',Param.preRowsTime,...
        %         'postRowsTime',Param.postRowsTime,...
        %         'prePhediTime',Param.prePhediTime,...
        %         'postPhediTime',Param.postPhediTime,...
        %         'scratchRegionMeters',Param.scratchRegionMeters,...
        %         'ExpDate',ExpDate,...
        %         'ExpHour',ExpHour...
        %         );
        %     %--- update the data struct
        %     exp_details=expDetailsRead(exp_dir);
        %     f = fieldnames(exp_details);
        %     for i = 1:length(f)
        %         ExperimentData.(f{i}) = exp_details.(f{i});
        %     end
        %     DataStruct.ExperimentData = ExperimentData;
        %
        %
        %     %% follow and analyze asperities
        %     disp('Definning and following asperities...');
        %     [vallyPairs,~,maxTrajectories,CM_Trajectories,skewnessMat,varMat,center_of_mass]...
        %         = Movie_followAsperities(RowOverTime);
        %
        %     %%
        %     fps = CamMetaStruct.FrameRate;        % frames/s
        %     postFrames = round(Param.postRowsTime*fps);
        %     preFrames = round(Param.preRowsTime*fps);
        %     %--- regulate number of post/pre frames
        %     if postFrames>CamMetaStruct.PostIms; postFrames=CamMetaStruct.PostIms; end
        %     if preFrames>(CamMetaStruct.NumIms-CamMetaStruct.PostIms); postFrames=CamMetaStruct.PostIms; end
        %
        %     spatialVec = (0:(size(RowOverTime,2)-1))/res;
        %
        %     %% save asperity data structure
        %     disp('saving asperity data...');
        %     analyzeAsperityStruct = struct(...
        %         'RowOverTime',RowOverTime,...
        %         'vallyPairs',vallyPairs,...
        %         'maxTrajectories',{maxTrajectories},...
        %         'CM_Trajectories',{CM_Trajectories},...
        %         'skewness',skewnessMat,...
        %         'var',varMat,...
        %         'center_of_mass',center_of_mass,...
        %         'spatialVec',spatialVec,...
        %         'CameraResolution',res,...
        %         'timeCount',timeCount,...
        %         'preFrames',preFrames,...
        %         'postFrames',postFrames...
        %         );
        %
        %     %--- update the data struct
        %     DataStruct.AsperityData = analyzeAsperityStruct;
        %
        %     %% get Big Picture data
        %     SubFolders = my_dir(exp_dir);
        %     if (sum(strcmpi(SubFolders,'Ph'))+sum(strcmpi(SubFolders,'PhBig')))==2
        %         BigPicRowOverTimeStruct = phantom_BigPicStruct(exp_dir, eventNum,Param.preRowsTime,1e-10,Param.postRowsTime,1,1,'all');
        %     else
        %         BigPicRowOverTimeStruct = IDT_PlotRowOverTime(exp_dir,eventNum,Param.preRowsTime,0,Param.postRowsTime,0,1:CamMetaStruct.ImageHeight,[],[],0);
        %     end
        %     %--- update the data struct
        %     DataStruct.BigPicRotStruct = BigPicRowOverTimeStruct;
                
        %% follow phedis
        
        %--- define frame to calc phedis on:
        prePhediFrames = round(Param.prePhediTime*CamMetaStruct.FrameRate);
        postPhediFrames = round(Param.postPhediTime*CamMetaStruct.FrameRate);
        firstFrame4Phedi = preFrames-prePhediFrames;
        lastFrame4Phedi = preFrames+postPhediFrames;
        if firstFrame4Phedi<1; firstFrame4Phedi=1; end;
        if lastFrame4Phedi>size(RowOverTime,1); lastFrame4Phedi=size(RowOverTime,1); end;
        
        %--- get phedis initial locations
%         FirstRowSig = RowOverTime(firstFrame4Phedi,:);
        FirstRowSig = RowOverTime_normalized(firstFrame4Phedi,:);
        phedisInitialLocsInPix = phedi_getPhediInitialLocs(FirstRowSig,Param.smooth4Phedi, Param.PhediSafetyDist, res, Param.scratchWidth);
        
        %--- calculate coarse shift
%         [coarseShiftsVector, shiftedRowOverTime] = Movie_findCoarseShiftsInRowOverTime(RowOverTime,center_of_mass);
        if Param.coarseShift
            [coarseShiftsVector, shiftedRowOverTime] = Movie_findCoarseShiftsInRowOverTime(RowOverTime_normalized,center_of_mass);
        else
            shiftedRowOverTime = RowOverTime_normalized;
            coarseShiftsVector = zeros(size(center_of_mass,1),1);
        end
        
        %--- follow phedis
        calShiftTimer = tic;
        if (~phantom_isPrecursor(DataStruct.BigPicRotStruct)) && Param.skipPrecursor
%             [frameCount, PhediLocationPixShifted, slopeIncline, StrechingFactorsMat] = Movie_phedi_calculate_continues_shift(...
%                 shiftedRowOverTime, firstFrame4Phedi, lastFrame4Phedi, phedisInitialLocsInPix,Param.numOOM,Param.phediRelative2first,Param.QuickAndDirty);
            

            %             [frameCount, PhediLocationPixShifted, slopeIncline, StrechingFactorsMat, SlopeRegions] = Movie_phedi_calculate_continues_shift(...
            %                 shiftedRowOverTime, firstFrame4Phedi, lastFrame4Phedi, phedisInitialLocsInPix,Param.numOOM,Param.phediRelative2first,slopePenalty, QuickAndDirty,[],interpP);
            [frameCount, PhediLocationPixShifted, slopeIncline, StrechingFactorsMat, SlopeRegions] = Movie_phedi_calculate_continues_shift_total(...
                shiftedRowOverTime, firstFrame4Phedi, lastFrame4Phedi, phedisInitialLocsInPix,...
                Param.followMethod,Param.numOOM, Param.slopePenalty, Param.expandLowBound, Param.expandHighBound, Param.phediRelative2first);
            
            relevantCoarseShifts = coarseShiftsVector(firstFrame4Phedi:lastFrame4Phedi);
            if nnz(coarseShiftsVector>0)
                PhediLocationPix = bsxfun(@plus,PhediLocationPixShifted,relevantCoarseShifts)-min(coarseShiftsVector(coarseShiftsVector>0));
            else
                PhediLocationPix = bsxfun(@plus,PhediLocationPixShifted,relevantCoarseShifts);
            end
            
            
            T_shiftTimer = toc(calShiftTimer);
            disp(['time for calculating phedi shifts = ',num2str(T_shiftTimer)]);
            
            %--- transform units: (minus 1 since the meters scale starts from 0 and the
            %pixel scale starts from 1)
            PhediLocation = (PhediLocationPix-1)/res;
%             measuredPhedisFromPlot = (phedisInitialLocsInPix-1)/res;
            measuredPhedisFromPlot = DataStruct.AsperityData.spatialVec(phedisInitialLocsInPix);
            timeVec= timeCount(frameCount);
            
            %% find phedis locaion
            [bigPicLocationPix, ~] = Movie_syncLocationOf2Cameras(BigPicRotStruct.t,BigPicRotStruct.DataMatNorm,...
                timeCount,RowOverTime,[],[],1);
%             PhotoLocation = Param.defPhotoLoc;
            PhotoLocation = BigPicRotStruct.x(bigPicLocationPix);
            
            [x_tips,t_tips] = phedi_find_t_tips_wPhotoLocation(BigPicRotStruct,PhotoLocation,DataStruct.AsperityData.spatialVec,...
                                    PhediLocation,measuredPhedisFromPlot,timeVec);
            
            
            %% phedi t_tips
            %--- find crack velocity in this region using big picture
            deltaX = 0.001; %region in which to integrate velocity (in meters)
            relevantRegionLogical = BigPicRotStruct.frontVelLoc_interpM>(Param.defPhotoLoc-deltaX) & BigPicRotStruct.frontVelLoc_interpM<(Param.defPhotoLoc+deltaX);
            Cf_MetPerSec = mean(BigPicRotStruct.frontVel_interpMperS(relevantRegionLogical));
            
            %--- get t_tips
%             CfPixPerFrame = Cf_MetPerSec*DataStruct.ExperimentData.res./DataStruct.CamMeta.FrameRate;
%             [intersectionFrame, tFrontPassing_Frms, crossLocInPix] = phedi_find_t_tips_wCf(...
%                 RowOverTime,CfPixPerFrame,phedisInitialLocsInPix,PhediLocation,firstFrame4Phedi:lastFrame4Phedi);
%             t_tips = interp1(1:size(RowOverTime,1),timeCount,tFrontPassing_Frms);
            CfCalcFrom = 'big pic photo lcoation';
            
            
            %--- get t_mins_t_tip
            t_mins_t_tip = bsxfun(@minus,timeVec(:),t_tips(:)');
            
            
            %% save Phedi data
            disp('saving phedi data...');
            analyzePhediStruct = struct(...
                'PhediLocation',PhediLocation,...
                'timeVec',timeVec,...
                't_mins_t_tip',t_mins_t_tip,...
                'measuredPhedisFromPlot',measuredPhedisFromPlot,...
                'slopeIncline',slopeIncline,...
                'firstFrame4Phedi',firstFrame4Phedi,...
                'lastFrame4Phedi',lastFrame4Phedi,...
                'phedisInitialLocsInPix',phedisInitialLocsInPix,...
                'PhediLocationPix',PhediLocationPix,...
                'StrechingFactorsMat',StrechingFactorsMat,...
                'SlopeRegions',SlopeRegions,...
                't_tips',t_tips,...
                'x_tips',x_tips,...
                'PhotoLocation',PhotoLocation,...
                'Cf',Cf_MetPerSec,...
                'CfCalcFrom',CfCalcFrom...
                );
            
            %--- update dataStruct
            DataStruct_WPhedis.PhediData = analyzePhediStruct;
        else
            DataStruct_WPhedis.PhediData = 'Precursor, phedi analysis skipped';
            disp('Precursor, phedi analysis skipped');
        end
        
        %% take time
        T_iteration = toc(phediAnalysisTimer);
        disp(['Elapsed time for this phedi analysis is ',num2str(T_iteration),' seconds']);
        fprintf('\n');
        fprintf('\n');
        
    catch
        %--- save the workspace in case of error
        w = whos;
        for a = 1:length(w)
            theWorkspace.(w(a).name) = eval(w(a).name);
        end
        DataStruct_WPhedis = struct;
        DataStruct_WPhedis.theWorkspace = theWorkspace;
        putvar('DataStruct');
        
    end
    
end