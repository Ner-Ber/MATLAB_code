function [analyzeDataStruct] = AnalyzeSingleMovie(figures_folder_path)
% 1. open movie in matlab movie player
% 2. press OK and enter event frame
% 3. open frame and mark asperities, save location data.

    %% parameters:
    bitDepth = 12;

    %% display movie and pick frame
    all_frames = createBWmovie(figures_folder_path,bitDepth);   % create the movie
    MovieDisplay = implay(all_frames);      % display movie so usr can find frame
    uiwait(MovieDisplay.Parent);            % wait for movie dialog to close in order to continue
    frameNum = str2double(cell2mat(inputdlg('enter frame of event')));      % ask and record frame number of event from usr.
    MovieDisplay = implay(all_frames);      % display movie so usr can find frame
    
    
    %% choose asperities areas
    coor_cell = {};
    measureFig = figure;
    imshow(all_frames(:,:,frameNum-2))
    h = imrect;
    
    %--- create buttons for UI
    uicontrol('Style', 'pushbutton', 'String', 'Done',...
            'Position', [20 20 100 40],...
            'Callback', 'close(gcf)');
        
    uicontrol('Style', 'pushbutton', 'String', 'start marking',...
        'Position', [130 20 100 40],...
        'Callback',@markAndMeasureRect);
   
    %--- terminate figure containing UI
    uiwait(gcf);         % wait for command and uiresume
    
    %--- convert coors measured into relevant units (integers and format of
    %[x_top y_top x_bottom y_bottom])
    coor_cell = convertCoors(coor_cell);
    
    
    %% create structure for export
    analyzeDataStruct = cell2struct({figures_folder_path, frameNum, coor_cell},...
        {'MoviePath','EventFrame','CoordinateCell'}, 2);
    
    
    
    
    
%% nested functions
%--- this function will execute measurment and record of rectangle:
    function markAndMeasureRect(src,event)
        while true
            p = wait(h);
            rectangle('Position', p, 'LineWidth',2, 'EdgeColor','r');
            coor_cell = cat(1, coor_cell, p);
        end
    end

%--- function that converts given rectangles into the required format
    function coor_cellConverted = convertCoors(coor_cell)
        coor_cellConverted = cell(size(coor_cell));
        for i = 1:length(coor_cell)
            c = coor_cell{i};
            coor_cellConverted{i} = [floor(c(1)), floor(c(2))+ceil(c(4)), floor(c(1))+ceil(c(3)), floor(c(2))];
        end
    end


end





