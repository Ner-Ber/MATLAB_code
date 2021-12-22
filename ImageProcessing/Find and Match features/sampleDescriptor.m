function desc = sampleDescriptor(im,pos,descRad)
% SAMPLEDESCRIPTOR Sample a MOPS?like descriptor at given positions in the image.
% Arguments:
% im ? nxm grayscale image to sample within.
% pos ? A Nx2 matrix of [x,y] descriptor positions in im.
% descRad ? ”Radius” of descriptors to compute (see below).
% Returns:
% desc ? A kxkxN 3?d matrix containing the ith descriptor
% at desc(:,:,i). The per?descriptor dimensions kxk are related to the
% descRad argument as follows k = 1+2*descRad.


    descNum = size(pos, 1);
    k = 1+2*descRad;
    desc = zeros(k, k, descNum);

    for i = 1:descNum
        x = pos(i,1)-descRad:pos(i,1)+descRad;
        y = (pos(i,2)-descRad:pos(i,2)+descRad)';
        [X,Y] = meshgrid(x,y);
        window = interp2(im,X,Y);           % currently accepts the 3rd level Gaussian pyramid!!!
        window(isnan(window)) = 0;
        Mean = mean(window(:));
        Norm = norm(window(:)-Mean);
        desc_current = (window - Mean)/Norm;
        desc_current(isnan(desc_current)) = 0;
        desc(:,:,i) = desc_current;
    end
    
end