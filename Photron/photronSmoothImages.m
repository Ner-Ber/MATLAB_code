function smoothedImage = photronSmoothImages(InitialImage, smt)
%smoothedImage = photronSmoothImages(InitialImage, smt)
% smoothedImage will smooth the InitialImage according to:
% if smt is a scalar, will average from all directions up to length of smt
% value.
% if smt is two element vector, will average in each direction according to
% legnth specified in vector [x_length y_length]

smt = fliplr(smt);
filter = ones(smt);
filter = filter/numel(filter);
smoothedImage = imfilter(InitialImage,filter,'replicate');