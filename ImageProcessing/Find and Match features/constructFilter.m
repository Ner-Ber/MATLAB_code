function filter = constructFilter(filterSize)
%constructFilter construct a 1D gaussian filter by a given size of
%filterSize

    a = [1 1];
    filter = a;
    for i = 1:(filterSize - 2)
        filter = conv(filter, a);
    end
    filter = filter/sum(filter);

end
