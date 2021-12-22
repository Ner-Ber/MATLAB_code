% props = expMeasureProperties(Exp)
%
% expMeasureProperties will measure properties that are relevant for all
% events in the experiment and save them to a file in the experiment's
% folder named 'exp_measures.txt'.
% In case such file already exists expMeasureProperties will read it and
% load the data.
%
% The option 'writeNew' will write a new such file even if one exists in
% the experiment's folder.

function expMeasures = expMeasureProperties(exp_dir,varargin)
    
    %% set defaults and parameters
    [writeNew, events,...
        pre_timeBig,    post_timeBig,   smtBig,     lineNumBig,...
        pre_timeSmall,  post_timeSmall,             lineNumSmall,   scratchDist...
        ] = setDefaults4function(varargin,...
        false, 'all',...
        2e-3,   4e-3,   0,  'all',...
        2e-3,   4e-3,       'all',      200e-6...
        );
    
    
    %% does the file already exist?
    Names = dir_names(exp_dir);
    existIndictr = sum(cellfun(@(A) strcmp(A,'exp_measures.txt'),Names));
    
    if existIndictr==0 || writeNew
        disp('creating a new ''exp_measures.txt'' file...');
        %% some parameters and data
        %-- events vector
        outdata = phantomReadMeta([exp_dir,'\Ph']);
        NumEvents = outdata.NumEvents;
        if strcmpi(events,'all')
            eventsVec = 1:NumEvents;
        else
            eventsVec = events(events<NumEvents & events>1);
        end
        
        %-- exp details
        exp_details = expDetailsRead(exp_dir);
        
        %% get big picture details
        blockEdgesCorrdinates = IDT_findEdgePixlesAllEvents(exp_dir,...
            events,pre_timeBig,[],post_timeBig,[],smtBig,lineNumBig);
        
        xx = linspace(0,exp_details.UpperBlockLength*1e-3,abs(diff(blockEdgesCorrdinates)));
        resBigPic = abs(xx(2)-xx(1)); % meters per pixel
        
        
        %% get small Pic details
        resSmall_Vec = [];
        preCours_Vec = [];
        PhotoLocation_Vec = [];
        for i = eventsVec
            try
                %-- small Pic ROT
                [RowOverTime_i, timeCount_i] = phantom_getRowOverTime_and_timeCount(exp_dir, i, pre_timeSmall, post_timeSmall, lineNumSmall);
                %-- Big Pic struct
                SubFolders = my_dir(exp_dir);
                
                if (sum(strcmpi(SubFolders,'Ph'))+sum(strcmpi(SubFolders,'PhBig')))==2
                    ims = phantomGetLines_BigPic(exp_dir,i,-pre_timeBig,1e-10,post_timeBig,1,1,'all',lineNumBig);
                else
                    if strcmpi(lineNumBig,'all')
                        CamMetaStruct = CameraMetaAllCams(exp_dir);
                        bigLinesIDT = 1:CamMetaStruct.ImageHeight;
                    else
                        bigLinesIDT = lineNumBig;
                    end
                    ims = IDT_imagesWithData(exp_dir,i,pre_timeBig,0,post_timeBig,0,bigLinesIDT);
                end
                RowOverTimeBig_i = ims.lines;
                %--- normalized ROT:
                DataMatNorm_i = bsxfun(@rdivide,RowOverTimeBig_i,RowOverTimeBig_i(1,:));
                alongFrameFixed = (1:size(RowOverTimeBig_i,2))-blockEdgesCorrdinates(1);
                xTemp = alongFrameFixed*resBigPic;
                
                
%                 if (sum(strcmpi(SubFolders,'Ph'))+sum(strcmpi(SubFolders,'PhBig')))==2
%                     %                     BigPicROT_i = phantom_BigPicStruct(exp_dir, i,pre_timeBig,1e-10,post_timeBig,1,1,'all');
%                     
%                     ims = phantomGetLines_BigPic(exp_dir,i,-pre_timeBig,1e-10,post_timeBig,1,1,'all',lineNumBig);
%                     RowOverTimeBig_i = ims.lines;
%                     %--- normalized ROT:
%                     DataMatNorm_i = bsxfun(@rdivide,RowOverTimeBig_i,RowOverTimeBig_i(1,:));
%                     alongFrameFixed = (1:size(RowOverTimeBig_i,2))-blockEdgesCorrdinates(1);
%                     xTemp = alongFrameFixed*resBigPic;
%                     
%                 else
%                     if strcmpi(lineNumBig,'all')
%                         CamMetaStruct = CameraMetaAllCams(exp_dir);
%                         bigLinesIDT = 1:CamMetaStruct.ImageHeight;
%                     else
%                         bigLinesIDT = lineNumBig;
%                     end
%                     ims = IDT_imagesWithData(exp_dir,i,pre_timeBig,0,post_timeBig,0,bigLinesIDT);
%                     RowOverTimeBig_i = ims.lines;
%                     %--- normalized ROT:
%                     DataMatNorm_i = bsxfun(@rdivide,RowOverTimeBig_i,RowOverTimeBig_i(1,:));
%                     alongFrameFixed = (1:size(RowOverTimeBig_i,2))-blockEdgesCorrdinates(1);
%                     
%                     CamMetaStruct = CameraMetaAllCams(exp_dir);
%                     BigPicROT_i = IDT_PlotRowOverTime(exp_dir,i,pre_timeBig,0,post_timeBig,0,bigLinesIDT,[],[],0);
%                 end
            catch
                continue
            end
            %-- get photoLocation
            [bigPicLocationPix, ~] = Movie_syncLocationOf2Cameras(...
                ims.t, DataMatNorm_i, timeCount_i, RowOverTime_i);
            if isnan(bigPicLocationPix)
                PhotoLocation_i = nan;
            else
                PhotoLocation_i = xTemp(bigPicLocationPix);
            end
            PhotoLocation_Vec = cat(1,PhotoLocation_Vec,PhotoLocation_i);
            %-- mark if precoursor for photolocation measurement
            preCours_i = phantom_isPrecursor(DataMatNorm_i);
            preCours_Vec = cat(1,preCours_Vec,preCours_i);
            %-- save small resolution
            resSmall_i = phantom_getResFromScrathces(RowOverTime_i, scratchDist);
            resSmall_Vec = cat(1,resSmall_Vec,resSmall_i);
        end
        
        %-- calculate small resolution
        resSmallPic = mean(resSmall_Vec);
        %-- calculate photolocation omiting nans and precoursors
        photoLocation = mean(PhotoLocation_Vec(~preCours_Vec),'omitnan');
        
        %% write to a file
        data2write = [min(blockEdgesCorrdinates), max(blockEdgesCorrdinates), resBigPic, resSmallPic, photoLocation];
        fileID = fopen([exp_dir,'\exp_measures.txt'],'w');
        fprintf(fileID,'%.14f\n',data2write);
        
        %-- add comments:
        comments = {'';'';'';'';'';'';'';...
            'comments:';...
            'row 1:';...
            'pixel in the big picture indicating left side of the upper block (computer side in system)';...
            'row 2:';...
            'pixel in the big picture indicating right side of the upper block (engine side in system)';...
            'row 3:';...
            'big picture resolution (m/pix)';...
            'row 4:';...
            'small picture resolution (pix/m)';...
            'row 5:';...
            'pixel where the center of the small camera frame is located on the big picture (m)';...
            '';'';'';...
            'created with the function ''expMeasureProperties'''};
        fprintf(fileID,'%s\n', comments{:});
        
        
        fclose(fileID);
        
    elseif existIndictr==1
        %% in case the file exists and needs to be read
        path_exp_measures=[exp_dir '\exp_measures.txt'];
        fileID = fopen(path_exp_measures);
        BigBlockLeft = str2double(fgetl(fileID));
        BigBlockRight = str2double(fgetl(fileID));
        resBigPic = str2double(fgetl(fileID));
        resSmallPic = str2double(fgetl(fileID));
        photoLocation = str2double(fgetl(fileID));
        blockEdgesCorrdinates = [BigBlockLeft, BigBlockRight];
        fclose(fileID);
    else
        %% in case somehow there's more than one file (?!)
        error(['more than 1 ''exp_measures.txt'' file in experiment ',exp_dir]);
    end
    
    %% write the data to a structure
    expMeasures = struct(...
        'blockEdgesCorrdinates',    blockEdgesCorrdinates,...
        'resBigPic',                resBigPic,...
        'resSmallPic',              resSmallPic,...
        'photoLocation',            photoLocation...
        );
end