function [pyr, filter] = LaplacianPyramid(im, maxLevels, filterSize)
%LaplacianPyramid constructs a Laplacian pyramid of a given image.
% im - a grayscale image with double values in [0,1] 
% maxLevels - the maximal number of levels in the resulting pyramid.
% filterSize - the size of the Gaussian filter (an odd scalar that
% represents a squared filter) to be used in constructing the pyramid filter.

    h = pyrHeight(im, maxLevels);
    %% Gaussian pyr
    Gpyr = GaussianPyramid(im, maxLevels, filterSize);
    
    %% construct filter
    filter = constructFilter(filterSize);
    
    %% construct pyramid
    pyr = cell(h,1);
    for i = 1:h-1
        S = size(Gpyr{i});
        pyr{i} = Gpyr{i} - expandIm(Gpyr{i+1}, filterSize, S(1), S(2));
    end
    pyr{h} = Gpyr{h};

end


