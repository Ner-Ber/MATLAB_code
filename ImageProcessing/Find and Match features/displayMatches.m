function displayMatches(im1,im2,pos1,pos2,inliers)
% DISPLAYMATCHES Display matched pt. pairs overlayed on given image pair.
% Arguments:
% im1,im2 ? two grayscale images
% pos1,pos2 ? Nx2 matrices containing n rows of [x,y] coordinates of matched
% points in im1 and im2 (i.e. the i’th match’s coordinate is
% pos1(i,:) in im1 and and pos2(i,:) in im2).
% inliers ? A kx1 vector of inlier matches (e.g. see output of
% ransacHomography.m)

    S = size(im1);
    IM = [im1; im2];
    pos2(:,2) = pos2(:,2)+S(1);

    outpos1 = pos1;
    outpos2 = pos2;
    outpos1(inliers,:) = [];
    outpos2(inliers,:) = [];
    pos1 = pos1(inliers,:);
    pos2 = pos2(inliers,:);

    % lines attachong inliers
    lineX = [pos1(:,1)'; pos2(:,1)'];
    lineY = [pos1(:,2)'; pos2(:,2)'];
    numPts = numel(lineX);
    lineX = [lineX; NaN(1,numPts/2)];
    lineY = [lineY; NaN(1,numPts/2)];

    % lines attachong outliers
    lineX_out = [outpos1(:,1)'; outpos2(:,1)'];
    lineY_out = [outpos1(:,2)'; outpos2(:,2)'];
    numPts_out = numel(lineX_out);
    lineX_out = [lineX_out; NaN(1,numPts_out/2)];
    lineY_out = [lineY_out; NaN(1,numPts_out/2)];

    % Display them!
    figure;
    imshow(IM);
    hold all;
    plot(lineX_out(:), lineY_out(:), 'b-');
    plot(outpos1(:,1),outpos1(:,2),'g.');
    plot(outpos2(:,1),outpos2(:,2),'g.');
    plot(lineX(:), lineY(:), 'y-');
    plot(pos1(:,1),pos1(:,2),'r.');
    plot(pos2(:,1),pos2(:,2),'r.');
    hold off;

end