function pos2 = applyHomography(pos1,H12)
% APPLYHOMOGRAPHY Transform coordinates pos1 to pos2 using homography H12.
% Arguments:
% pos1 ? An Nx2 matrix of [x,y] point coordinates per row.
% H12 ? A 3x3 homography matrix.
% Returns:
% pos2 ? An Nx2 matrix of [x,y] point coordinates per row obtained from
% transforming pos1 using H12.

    S = size(pos1);
    coors1 = [permute(pos1,[2 1]); ones(1,S(1))];
    coors2 = H12*coors1;
    pos2 = coors2(1:2,:)'./[coors2(3,:)',coors2(3,:)'];



end