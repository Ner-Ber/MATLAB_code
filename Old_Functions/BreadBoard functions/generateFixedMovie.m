function [fixed_frames, fixed_frames_RGB, all_frames] = generateFixedMovie(folderPath, movie_name, movieType,fix_blinking)
%generateFixedMovie will generate fixed frames from given blinking frames 
% folderPath - folder path of movie frames. if insereted empty string ('')
% a UI wondow will appear to ask for a folder path.
% movie_name - name of movie to be saved (can include path, USE SUFFIX!).
% is entered empty the movie will not be saved.
% movieType - string indicating which type of movie to save. use 'RGB' or
% 'BW' for each type. if inserted empty  ('') or other no movie will be
% saved.

%% create matrix of frames
    %---parameters
    bitDepth = 12;
    ColorMap = colormap(jet);
    crop = 1:8;

    switch folderPath
        case ''
            [all_frames, all_frames_RGB] = CreateMovieFrames(bitDepth, ColorMap, crop);
        otherwise
            [all_frames, all_frames_RGB] = CreateMovieFrames(bitDepth, ColorMap, crop, folderPath);
    end

%% fix blinking
    switch fix_blinking
        case 1
            %---parameters:
            FirstFrames = 30;
            intensityThreshold = 0.1;

            [fixed_frames, ~] = fixBlinking(all_frames, FirstFrames, intensityThreshold);
            [fixed_frames_RGB, ~] = fixBlinking(all_frames_RGB, FirstFrames, intensityThreshold);
   
        case 0
            fixed_frames = all_frames;
            fixed_frames_RGB = all_frames_RGB;
    end    


%% create movie

%----fix all the pixels that are more then 1 so the movie will work


%---switch to saving movie type
    switch movieType
        case 'BW'   %save the BW movie
            %----fix all the pixels that are more then 1 or les  than 0 so the movie will work
            movie_frames=fixed_frames;
%             d=logical(movie_frames>1);                %%find the pixels over 1
%             e=logical(movie_frames<0);                %%find the pixels under 0
%             movie_frames(d)=floor(movie_frames(d)); %%make them 1 for the movie
%             movie_frames(e)=ceil(movie_frames(e)); %%make them 0 for the movie
            saving_frames =  movie_frames;
        case 'RGB'	%save the RGB movie
            saving_frames = fixed_frames_RGB;
        otherwise   %no proper indication of movie type. no movie will be saved.
            movie_name = '';
%             warning('movieType not defined well. No movie will be saved');
    end
    
        
    %---save movie
    switch movie_name
        case ''
            return;
        otherwise
            videoMakerHandle = VideoWriter(movie_name);
            open(videoMakerHandle);
            for i = 1:size(fixed_frames,3)
                writeVideo(videoMakerHandle, saving_frames(:,:,i));
            end
            close(videoMakerHandle)
    end
end