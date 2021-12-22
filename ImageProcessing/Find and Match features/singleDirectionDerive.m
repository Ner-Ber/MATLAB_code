function Ii = singleDirectionDerive( im , dir )
%singleDirectionDerive derives a given image in a steted direction
% Arguments:
% im ? nxm grayscale image to derive
% dir - directoin of derivation. Consistant with matlab directions and axis
% dir='y' - row direction, y
% dir='x' - column direction, x


    ker = [1 0 -1];
    if strcmp(dir,'y')
        rot = 1;
    elseif strcmp(dir,'x')
        rot = 0;
    else
        error('Invalid "dir" input in singleDirectionDerive');
    end

    if rot
        ker = ker.';
    end

    Ii = conv2(im, ker, 'same');



end