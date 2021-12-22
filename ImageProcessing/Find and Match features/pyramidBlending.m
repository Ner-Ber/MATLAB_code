function imBlend = pyramidBlending(im1, im2, mask, maxLevels, filterSizeIm, filterSizeMask)
%#pyramidBlending blends two images according to a given mask
%#im1, im2 - are two input grayscale images to be blended.
%#mask - is a binary mask containing 1’s and 0’s representing which parts of im1 and im2 should appear in
% the resulting imBlend.
%#maxLevels - is the maxLevels parameter you should use when generating the Gaussian and Laplacian
% pyramids.
%#filterSizeIm - is the size of the Gaussian filter (an odd scalar that represents a squared filter) which
% defining the filter used in the construction of the Laplacian pyramids of im1 and im2.
%#filterSizeMask - is the size of the Gaussian filter(an odd scalar that represents a squared filter) which
% 4 defining the filter used in the construction of the Gaussian pyramid of mask.
%#Note that im1, im2 and mask should all have the same dimensions.


    [L1, filter] = LaplacianPyramid(im1, maxLevels, filterSizeIm);
    [L2, ~] = LaplacianPyramid(im2, maxLevels, filterSizeIm);
    Gm = GaussianPyramid(mask, maxLevels, filterSizeMask);
    h = length(L1);
    Lout = cell(h,1);
    
    for k = 1:h 
        Lout{k} = Gm{k}.*L1{k} + (ones(size(Gm{k}))-Gm{k}).*L2{k};
    end

    imBlend = LaplacianToImage(Lout, filter, ones(h,1));


end



function img = LaplacianToImage(lpyr, filter, coeffMultVec)
%LaplacianToImage reconstructs the original image from a laplacian pyramid

    h = length(lpyr);
    filterSize = length(filter);
    mul_pyr = cell(h,1);
    for i = 1:h
        mul_pyr{i} = lpyr{i}*coeffMultVec(i);
    end

    %% construct fullSizePyr
    recImPyr = cell(h,1);
    recImPyr{h} = mul_pyr{h};
    for i = h-1:-1:1
        S = size(mul_pyr{i});
        recImPyr{i} = mul_pyr{i} + expandIm(recImPyr{i+1}, filterSize, S(1), S(2));
    end

    img = recImPyr{1};

end