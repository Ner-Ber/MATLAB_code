function [H12,inliers] = ransacHomography(pos1,pos2,numIters,inlierTol)
% RANSACHOMOGRAPHY Fit homography to maximal inliers given point matches
% using the RANSAC algorithm.
% Arguments:
% pos1,pos2 ? Two Nx2 matrices containing n rows of [x,y] coordinates of
% matched points.
% numIters ? Number of RANSAC iterations to perform.
% inlierTol ? inlier tolerance threshold.
% Returns:
% H12 ? A 3x3 normalized homography matrix.
% inliers ? A kx1 vector where k is the number of inliers, containing the indices in pos1/pos2 of the maximal set of
% inlier matches found.

    N = size(pos1, 1);
    E_rec = [-inf, 0, 0, 0, 0];

    for i = 1:numIters
        J = randperm(N, 4);
        H12 = leastSquaresHomography(pos1(J,:),pos2(J,:));
        if isempty(H12)
            continue
        end

        pos1_t = applyHomography(pos1,H12);
        Ej = sum(abs(pos1_t - pos2).^2,2);
        E_tot = length(find(Ej<inlierTol));
        if E_tot > E_rec(1)
            E_rec = [E_tot,J];
        end   
    end
    Jin = E_rec(2:end);
    H12 = leastSquaresHomography(pos1(Jin,:),pos2(Jin,:));
    pos1_t = applyHomography(pos1,H12);
    Ej = sum(abs(pos1_t - pos2).^2,2);
    inliers = find(Ej<inlierTol);
    H12 = leastSquaresHomography(pos1(inliers,:),pos2(inliers,:));

end