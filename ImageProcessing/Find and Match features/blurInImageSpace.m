function blurImage = blurInImageSpace(inImage,kernelSize)
%blurInImageSpace function that performs image blurring using 2D
%convolution between the image f and a gaussian kernel g
% inImage - is the input image to be blurred (grayscale double image).
% kernelSize - is the size of the gaussian in each dimension (one odd integer).
% blurImage - is the output blurry image (grayscale double image).


    g = createGaussianAprox(kernelSize);        % creating the kernel g
    blurImage = conv2(inImage, g, 'same');      % using convolution to blur


end

