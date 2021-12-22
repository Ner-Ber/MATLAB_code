function [ imQuant, error] = quantizeImage(imOrig, nQuant, nIter)

%quantizeImage   performs optimal quantization of a given grayscale or RGB image. 
% imOrig - is the input grayscale or RGB image to be quantized.
% nQuant - is the number of intensities your output imQuant image should have.
% nIter - is the maximum number of iterations of the optimization procedure. (May converge earlier.)


    % work on proper layer depending on RGB/Grayscale
    if size(imOrig, 3) == 3
        imOrigYIQ = rgb2ntsc(imOrig);
        im = imOrigYIQ(:,:,1);
    else
        im = imOrig;
    end

    im = uint8(im*255);
    HistOrig = imhist(im);
    CumHist = cumsum(HistOrig);
    pix_dens = numel(im)/nQuant;                            % define number of pixels per quanta (rough approx)
    z = cumsum(accumarray(floor(CumHist/pix_dens) + 1,1));  % create a initial vector of Quanta's separators
    z = [1;z(1:(end-1))];                                   % define boundary condition #1 beggining
    z(end) = 256;                                           % define boundary condition #2 end
    q = [];
    Er_iSq = [];
    error = [];
    z_former = -1.*ones(length(z),1);                       % defining former sparation vec just for the prot
    j = 1;
    while j <= nIter && ~isequal(z_former,z)
        for i = 1:nQuant
            % define the optimal grayscale-tone according to current z vec
            q(i) = sum(HistOrig(z(i):z(i+1))'.*(z(i):z(i+1)))/(sum(HistOrig(z(i):z(i+1))));
            q(i) = round(q(i));
            Er_iSq(i) = sum((((z(i):z(i+1))-q(i)).^2)'.*HistOrig(z(i):z(i+1)));
        end
        error(j) = sum(Er_iSq);
        z_former = z;    
        for i = 2:nQuant
            z(i) = (q(i)+q(i-1))/2; 
        end
        z = round(z);    
        j = j+1;         
    end
    error = error';
   
    % creating the lookup table
    LUT = zeros(length(HistOrig),1);
    for i = 1:nQuant;
        LUT(z(i):end) = q(i);
    end
    
    % applying the LUT
    imQuant = double(LUT(im+1)/255);
    if size(imOrig, 3) == 3
        imQuant(:,:,2:3) = imOrigYIQ(:,:,2:3);
        imQuant = ntsc2rgb(imQuant);
    end
    
    % diplay
%     figure('Name', 'quantization error along iterations'); plot(error);
%     figure('Name', 'Original Image'); imshow(imOrig);
%     figure('Name', 'Quantizesed Image'); imshow(imQuant);

end

