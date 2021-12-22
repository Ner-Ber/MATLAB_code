function panorama = renderPanorama(im,H,Xmin,Ymin,Ipano,masksCell)
% RENDERPANORAMA Renders a set of images into a combined panorama image.
% Arguments:
% im ? Cell array of n grayscale images.
% H ? Cell array of n 3x3 homography matrices transforming the ith image
% coordinates to the panorama image coordinates.
% Returns:
% panorama ? A grayscale panorama image composed of n vertical strips that
% were backwarped each from the relevant frame im{i} using homography H{i}.

    %% create a canvas and compute strips
    M = length(im);
    xmins = zeros(1,M);
    xmaxs = zeros(1,M);
    ymins = zeros(1,M);
    ymaxs = zeros(1,M);
    HomIm = cell(M,1);

    %% compute warped images and create masks
    for i = 1:M
        [HomIm{i}, x2, y2] = warpImage(im{i}, H{i});
        xmins(i) = min(x2(:));
        xmaxs(i) = max(x2(:));
        ymins(i) = min(y2(:));
        ymaxs(i) = max(y2(:));
    end

    %% creating same-sized images
    for i = 1:M
        %% padding images 
        Padded1 = Ipano;
        Padded1((ymins(i)-Ymin+1):(ymaxs(i)-Ymin+1),...
               (xmins(i)-Xmin+1):((xmaxs(i)-Xmin+1))) = HomIm{i};
        HomIm{i} = Padded1;
    end

    maxLevels = 7;
    filterSizeIm = 25;
    filterSizeMask = 25;
    S = size(Ipano);
    im1 = [];

    %% blend
    for i = 1:M-1
        if isempty(im1)
            im1 = HomIm{i};
        end
        im2 = HomIm{i+1};
        imBlend = pyramidBlending(im1, im2,...
            masksCell{i}, maxLevels, filterSizeIm, filterSizeMask);

        im1 = imBlend;
    end

    im1 = im1(1:S(1), 1:S(2));  % make sure the panorama to the original size
    panorama = im1;

end