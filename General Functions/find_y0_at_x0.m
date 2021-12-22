function y0 = find_y0_at_x0(x0,y,x)
%y0 = FIND_Y0_AT_X0(x0,y,x)
%
%FIND_Y0_AT_X0 will find the value of the function y=f(x) at the point x0.

[~,II] = min(abs(x-x0));
y0 = y(II);

end