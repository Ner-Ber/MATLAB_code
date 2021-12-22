function pyr = GaussianPyramid(im, maxLevels, filterSize)
%GaussianPyramid constructs a Gaussian pyramid of a given image.
% im - a grayscale image with double values in [0,1] 
% maxLevels - the maximal number of levels in the resulting pyramid.
% filterSize - the size of the Gaussian filter (an odd scalar that
% represents a squared filter) to be used in constructing the pyramid filter.

    %% calculate hieght of pyramid
    h = pyrHeight(im, maxLevels);
    
    %% construct pyramid
    pyr = cell(h,1);
    pyr{1} = im;
    for i = 2:h
        blur_im = blurInImageSpace(pyr{i-1},filterSize);
        reduced_im = blur_im(1:2:end, 1:2:end);
        pyr{i} = reduced_im;
    end

end