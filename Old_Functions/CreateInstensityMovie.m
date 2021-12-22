function [all_frames, all_frames_RGB, videoMakerHandle]...
    = CreateInstensityMovie(MovieName, bitDepth, ColorMap, crop)
%plotPixelRowOverTime will plot the intensity change of a certain row in
%the images inserted as function of time.
% MovieName - name the movie and location to save. include suffix.
% bitDepth - bit depth of the movie frames
% ColorMap - colormap you'd like to use when displaying intensities
% crop - defines rows which will be displayed. insert 0 for no croppping.
% other wise insert the row number ot numbers you'd like to view. 

%---load images
[image_names,figures_folder_path,~] = uigetfile({'*.jpg;*.tif;*.png;*.gif;*.ibw','All Image Files';...
'*.*','All Files' },'select images to find edges (**same batch**)','MultiSelect','on');

%---create video writer
videoMakerHandle = VideoWriter(MovieName);
open(videoMakerHandle);

%---define crop condition
if isscalar(crop)
    if crop == 0
        cropBorder = size(imread(fullfile(figures_folder_path,image_names{1})),1);
    else
        cropBorder = crop;
    end
else
    cropBorder = crop;
end

%---read images and create write video
all_frames_RGB = [];
all_frames = [];
for FramIdx = 1:2:length(image_names)
    currentImage = double(imread(fullfile(...
        figures_folder_path,image_names{FramIdx})))./2^(bitDepth);
    currentFrame = currentImage(1:cropBorder,:);
    RGB_frame = ind2rgb(gray2ind(currentFrame),ColorMap);
    all_frames_RGB = cat(4, all_frames_RGB, RGB_frame);
    all_frames = cat(3, all_frames, currentFrame);
    writeVideo(videoMakerHandle, RGB_frame);
end
close(videoMakerHandle);

end