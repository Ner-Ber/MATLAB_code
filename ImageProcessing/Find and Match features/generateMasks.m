function [masksCell,Xmin,Ymin,Ipano] = generateMasks(im,H)
% generateMasks creates a cell array of masks that fits a set of RGB images
% intended to be stiched
% the function retrieves more than masks:
% Xmin/Ymin - the top coenrs in the panorama coordinate system
% Ipano - an empty ("0") canvas the size of the intended panorama

    %% create a canvas and compute strips
    M = length(im);
    xmins = zeros(1,M);
    xmaxs = zeros(1,M);
    ymins = zeros(1,M);
    ymaxs = zeros(1,M);
    HomIm = cell(M,1);
    masksCell = cell(M-1,1);
    
    %% compute warped images and create masks
    for i = 1:M
        [HomIm{i}, x2, y2] = warpImage(im{i}, H{i});
        xmins(i) = min(x2(:));
        xmaxs(i) = max(x2(:));
        ymins(i) = min(y2(:));
        ymaxs(i) = max(y2(:));
    end

    Xmin = min(xmins);
    Ymin = min(ymins);
    Xmax = max(xmaxs);
    Ymax = max(ymaxs);
    pcenters = mean([xmins; xmaxs]);
    dx = abs(Xmin - Xmax)+1;
    dy = abs(Ymin - Ymax)+1;
    P = round(pcenters-Xmin);
    Ipano = zeros(dy, dx);      % create a canvas

    %% creating same-sized images and pyramid blending
    for i = 1:M
        %% padding images and masks to fit the full-size panorama
        Padded1 = Ipano;
        Padded1((ymins(i)-Ymin+1):(ymaxs(i)-Ymin+1),...
               (xmins(i)-Xmin+1):((xmaxs(i)-Xmin+1))) = HomIm{i};

        HomIm{i} = Padded1;
    end
    
    %% create masks
    for i = 1:M-1
        im1 = HomIm{i};
        im2 = HomIm{i+1};
        %% create a strip for dynamic programming
        strip1 = im1(:,P(i):P(i+1));
        strip2 = im2(:,P(i):P(i+1));
        Mask_strip = dynProgram(strip1, strip2);
        Mask = padarray(Mask_strip, [0 P(i)-1 0], 1, 'pre');
        masksCell{i} = padarray(Mask, [0 dx-P(i+1) 0], 0, 'post');


    end

end