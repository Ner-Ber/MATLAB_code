function pos = HarrisCornerDetector( im )
% HARRISCORNERDETECTOR Extract key points from the image.
% Arguments:
% im - nxm grayscale image to find key points inside.
% pos - A nx2 matrix of [x,y] key points positions in im

    %% define variables
    kernelSize = 3;
    k = 0.04;

    %% Define M matrix elements
    Ix = singleDirectionDerive( im , 'x' );
    Iy = singleDirectionDerive( im , 'y' );
    Ix_sq = blurInImageSpace(Ix.^2,kernelSize);
    Iy_sq = blurInImageSpace(Iy.^2,kernelSize);
    IxIy = blurInImageSpace(Ix.*Iy,kernelSize);
    IyIx = blurInImageSpace(Iy.*Ix,kernelSize);

    %% clculate R = response image & give maxima coordinates
    DetM = Ix_sq.*Iy_sq - IxIy.*IyIx;
    TraceM = Ix_sq+Iy_sq;
    response = DetM -k*(TraceM.^2);
    result = nonMaximumSuppression( response );
    [row,col] = find(result);
    
    pos = [col, row];

end