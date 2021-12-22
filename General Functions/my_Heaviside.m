function y = my_Heaviside(x)
    Logic = x>=0;
    y = zeros(size(x));
    y(Logic) = 1;
end