function [x2, y2] = generateImageGrids(im1, H12)

% generateImageGrids genrates the grids needed for the interp2 in order to
% apply an homography transformation on a given image im1.
% syntax:
% [x2, y2] = generateImageGrids(im1, H12)
% 
% im1 - a given grayscale double image
% H12 - a 3x3 linear transformation
% x2, y2 - new coordiante sets

    S = size(im1);
    [x1,y1] = meshgrid(1:S(2),1:S(1));
    x1 = x1(:);
    y1 = y1(:);
    pos1 = [x1, y1];
    pos2 = applyHomography(pos1,H12);
    x2 = reshape(pos2(:,1),S);
    y2 = reshape(pos2(:,2),S);


end