function im = imReadAndConvert( filename, represantation )
%imReadAndConvert reads a given image file and converts it into a given representation.
%   filename - should be tha path of the file or it's name in the opened
%   folder. 
%   represantation - either 1 or 2 defining if the output should be either a grayscale
%   image (1) or an RGB image (2)

    
    inf = imfinfo(filename);                    % creating variable of inf for future use
    import = imread(filename);
    import = double(import)/255;
    if represantation == 1                      % converting to grayscale if needed
        if ~strcmp(inf.ColorType,'grayscale')   % using the info about the image
            im = rgb2gray(import);
        else
            im = import;
        end
    elseif represantation ==2
        im = import;
    end   


end

