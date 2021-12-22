function Htot = accumulateHomographies(Hpair,m)
% ACCUMULATEHOMOGRAPHY Accumulate homography matrix sequence.
% Arguments:
% Hpair ? Cell array of M?1 3x3 homography matrices where Hpair{i} is a
% homography that transforms between coordinate systems i and i+1.
% m ? Index of coordinate system we would like to accumulate the
% given homographies towards (see details below).
% Returns:
% Htot ? Cell array of M 3x3 homography matrices where Htot{i} transforms
% coordinate system i to the coordinate system having the index m.
% Note:
% In this exercise homography matrices should always maintain
% the property that H(3,3)==1. This should be done by normalizing them as
% follows before using them to perform transformations H = H/H(3,3).

    Hpair = cellfun(@(x) x./x(3,3), Hpair, 'UniformOutput', 0);
    M = length(Hpair)+1;
    Htot = cell(M,1);
    pre = m-1:-1:1;
    post = m:1:M-1;

    H = [];
    for i = pre
        if isempty(H);
            H = Hpair{i};
        else
            H = H*Hpair{i};
        end
        Htot{i} = H;
    end

    H = [];
    for i = post
        if isempty(H);
            H = inv(Hpair{i});
        else
            H = H/Hpair{i};
        end
        Htot{i+1} = H;
    end

    Htot{m} = eye(3);

end