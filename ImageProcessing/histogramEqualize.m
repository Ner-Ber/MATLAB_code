function [ imEq, HistOrig, histEq ] = histogramEqualize( imOrig )
% performs histogram equalization of a given grayscale or RGB image
%   imOrig - is the input grayscale or RGB double image with values in [0, 1].
%   imEq - is the equalized image. grayscale or RGB double image with values in [0, 1].
%    histOrig - is a 256 bin histogram of the original image (256x1 vector).
%    histEq - is a 256 bin histogram of the equalized image (256x1 vector)


    % work on proper layer depending on RGB/Grayscal
    if size(imOrig, 3) == 3
        imOrigYIQ = rgb2ntsc(imOrig);
        im = imOrigYIQ(:,:,1);
    else
        im = imOrig;
    end
    
    HistOrig = imhist(im);
    im = im*255;
    
    % creating the Lookup table
    CumHistOrig = cumsum(HistOrig);
    CumHistOrig_norm = CumHistOrig./(numel(im(:,:,1)));
    CumHistOrig_norm = CumHistOrig_norm*255;
    CumHistOrig_R = round(CumHistOrig_norm);
    
    % applying the LUT
    imEq = CumHistOrig_R(uint8(im+1));
    imEq = imEq/255;
    histEq = imhist(imEq);
    
    % converting back to RGB if needed
    if size(imOrig, 3) == 3
        imEq(:,:,2:3) = imOrigYIQ(:,:,2:3);
        imEq = ntsc2rgb(imEq);
    end
    
    % displaying
    figure('Name', 'Original Image'); imshow(imOrig);
    figure('Name', 'Equalized Image'); imshow(imEq);
end
