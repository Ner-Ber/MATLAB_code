function [pos,desc] = findFeatures(pyr)
% FINDFEATURES Detect feature points in pyramid and sample their descriptors.
% This function should call the functions spreadOutCorners for getting the keypoints, and
% sampleDescriptor for sampling a descriptor for each keypoint
% Arguments:
% pyr ? Gaussian pyramid of a grayscale image having 3 levels.
% Returns:
% pos ? An Nx2 matrix of [x,y] feature positions per row found in pyr. These
% coordinates are provided at the pyramid level pyr{1}.
% desc ? A kxkxN feature descriptor matrix.

    %% parameters
    n = 7;
    m = 7;
    radius = 2;
    descRad = 3;
    P = length(pyr);
    pos = [];
    desc = [];

    %% get descs from all layers of pyramid
    for i = 1:P
        pos_i = spreadOutCorners(pyr{i},n,m,radius);
        coor_P = pyrLevelCoorTrans(pos_i, i, P);
        coor_1 = pyrLevelCoorTrans(pos_i, i, 1);
        desc = cat(3, desc, sampleDescriptor(pyr{P},coor_P,descRad));
        pos = cat(1, pos, coor_1);
        
        
    end
end