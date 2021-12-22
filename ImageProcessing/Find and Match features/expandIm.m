function exp_im = expandIm(im, filterSize, m, n)
%expandIm expand a given image "im" by factor of n, by zero padding in
%between line and blurring by a given image size

% this function has been changed since the last ex3 so it will be able to
% expand in non 2^n resolutions

    canv = zeros(m,n);
    canv(1:2:m, 1:2:n) = im(1:length(1:2:m),1:length(1:2:n));
    filter = constructFilter(filterSize);
    blur_x = conv2(canv, 2*filter, 'same');
    exp_im = conv2(blur_x, 2*filter', 'same');

end