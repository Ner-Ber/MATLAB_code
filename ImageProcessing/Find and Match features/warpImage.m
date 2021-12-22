function [im_rot, x2, y2] = warpImage(im1, H12)
% warpImage warps imahes using a linear homography represented by the 3x3
% matrix H12,
% returns also a mesh grid in the new coordinate system

    [x2, y2] = generateImageGrids(im1, H12);
    xmin = floor(min(x2(:)));
    ymin = floor(min(y2(:)));
    xmax = ceil(max(x2(:)));
    ymax = ceil(max(y2(:)));
    [x2,y2] = meshgrid(xmin:xmax,ymin:ymax);

    pos_back = applyHomography([x2(:),y2(:)], inv(H12));
    X = reshape(pos_back(:,1), size(x2,1),size(x2,2));
    Y = reshape(pos_back(:,2), size(y2,1),size(y2,2));

    S = size(im1);
    [x1,y1] = meshgrid(1:S(2),1:S(1));

    im_rot = interp2(x1,y1,im1,X,Y,'linear');
    im_rot(isnan(im_rot)) = 0;

end