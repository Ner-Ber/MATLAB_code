function [DataStruct_WPhedis] = Movie_phedi_from_folder_2_data_parameterStruct(DataStruct)
    % [DataStruct] = Movie_phedi_from_folder_2_data_parameterStruct(experNum, eventNum,varargin)
    %
    % Movie_phedi_from_folder_2_data will read the images from the movie, then
    % claculate phedis and follow them end output cell arrays containig the
    % relevant data.
    %
    %
    %
    %   09-12-2018
    %       *) 'prePhediTime' and 'postPhediTime' are now compared to
    %       mean(t_tips) and not to the trigger time. 
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
        
        %% get photo location and time
        expMeasures = expMeasureProperties(DataStruct.ExperimentData.ExpHour);
        PhotoLocation = expMeasures.photoLocation;
        timeVec4space = (BigPicRotStruct.frontTime_interp)/BigPicRotStruct.fps;
        [~,bigPicLocationPix] = min(abs(BigPicRotStruct.x-PhotoLocation));
        centerTimeSec = timeVec4space(bigPicLocationPix);
        centerTimeFrms = round(centerTimeSec*fps);
            
%         [bigPicLocationPix, ~] = Movie_syncLocationOf2Cameras(BigPicRotStruct.t,BigPicRotStruct.DataMatNorm,...
%             timeCount,RowOverTime,[],[],1);
%         %             PhotoLocation = Param.defPhotoLoc;
%         if ~isnan(bigPicLocationPix)
%             PhotoLocation = BigPicRotStruct.x(bigPicLocationPix);
%             %-- get photo time in frames of the small pic for phedi analysis
%             timeVec4space = (BigPicRotStruct.frontTime_interp)/BigPicRotStruct.fps;
%             centerTimeSec = timeVec4space(bigPicLocationPix);
%             centerTimeFrms = round(centerTimeSec*fps);
%         else
%             PhotoLocation = Param.defPhotoLoc;
%             centerTimeFrms = 0;
%         end
        
        
        
        
        %% follow phedis
        
        %--- define frame to calc phedis on:
        prePhediFrames = round(Param.prePhediTime*CamMetaStruct.FrameRate)-centerTimeFrms;
        postPhediFrames = round(Param.postPhediTime*CamMetaStruct.FrameRate)+centerTimeFrms;
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
                Param.followMethod,Param.numOOM, Param.slopePenalty, Param.expandLowBound, Param.expandHighBound, Param.phediRelative2first, Param.signalSmooth4PhediTrack);
            
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
            
            %% find phedis location and time          
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