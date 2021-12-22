function all_frames = createBWmovie(ImagesPath,bitDepth)
%createBWmovie creates a simple black and white movie from the frames saved
%as images in the folder entered. 


    %---load images
    listing = dir([ImagesPath,'\*.tif']);
    image_names = {listing.name}';
    
    %---read images and create write video
    all_frames = [];
    for FramIdx = 1:length(image_names)
        currentImage = double(imread(fullfile(...
            ImagesPath,image_names{FramIdx})))./2^(bitDepth);
        all_frames = cat(3, all_frames, currentImage);
    end
end