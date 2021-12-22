function [all_frames, all_frames_RGB] = CreateMovieFrames(bitDepth, ColorMap, crop, figures_folder_path)
%CreateMovieFrames will create two matrices of movie frames in gray shades
%and in RGB.
% all_frames - black and white frames organized in a 3D matrix. dimensions
% 1&2 are y and x axis of the frame, dimension 3 is time progress.
% all_frames_RGB - RGB frames organized in a 4D matrix. dimensions
% 1&2 are y and x axis of the frame, dimension 3 is RGB layers (n this order)
% dimension 4 is time progress.
% MovieName - name the movie and location to save. include suffix.
% bitDepth - bit depth of the movie frames
% ColorMap - colormap you'd like to use when displaying intensities
% crop - defines rows which will be displayed. insert 0 for no croppping.
% other wise insert the row number ot numbers you'd like to view. 

    %---load images
    switch nargin
        case 3
        [image_names,figures_folder_path,~] = uigetfile({'*.jpg;*.tif;*.png;*.gif;*.ibw','All Image Files';...
        '*.*','All Files' },'select images to find edges (**same batch**)','MultiSelect','on');
        case 4
        listing = dir([figures_folder_path,'\*.tif']);
        image_names = {listing.name}';
    end


    %---define crop condition
    if isscalar(crop)
        if crop == 0
            cropBorder = 1:size(imread(fullfile(figures_folder_path,image_names{1})),1);
        else
            cropBorder = crop;
        end
    else
        cropBorder = crop;
    end

    %---read images and create write video
    all_frames_RGB = [];
    all_frames = [];
    for FramIdx = 1:length(image_names)
        currentImage = double(imread(fullfile(...
            figures_folder_path,image_names{FramIdx})))./2^(bitDepth);
        currentFrame = currentImage(cropBorder,:);
        RGB_frame = ind2rgb(gray2ind(currentFrame),ColorMap);
        all_frames_RGB = cat(4, all_frames_RGB, RGB_frame);
        all_frames = cat(3, all_frames, currentFrame);
    end

end