function x0 = find_x0_at_y0(y0,y,x,varargin)
%x0 = FIND_X0_AT_Y0(y0,y,x,pref)
%
%FIND_X0_AT_Y0 will find a location where the function y=f(x) gets the
%value of y0. this will be done by linear interpolation.
%
%INPUTS:
% y0 - the value of the function of where you want to find the location.
% y - vector defining the values of the function
% x - vector defining the location of the function
% pref - prefaration in case multiple values occur. this variable can get
% the values:
%	+) 'mean' - will return the mean value of all results (default)
%	+) 'left' - will return the value corrosponding to the smallest x val
%	+) 'right' - will return the value corrosponding to the largest x val
%   +) 'all' - returns a vector containing all crossing points

[pref] = setDefaults4function(varargin,'mean');

zeroed_y = y-y0;
%--- find where exact value meets:
exact_x0 = x(zeroed_y==0);
%--- find zero crossings:
positive = zeroed_y>0;
positive = positive(:)';
zeroCross = sum([positive(1:end-1);positive(2:end)],1)==1;  %++ here the index i represents weather there is crossing between indexes i to i+1
x1 = x(zeroCross);
x2 = x(logical([0 zeroCross(1:end-1)]));
y1 = zeroed_y(zeroCross);
y2 = zeroed_y(logical([0 zeroCross(1:end-1)]));
m = (y2(:)-y1(:))./(x2(:)-x1(:));
x_crossing = x1(:)-(1./m).*y1;

x0_vec = unique([exact_x0(:);x_crossing(:)]);
if strcmp(pref,'mean')
    x0 = mean(x0_vec);
elseif strcmp(pref,'left')
    x0 = min(x0_vec);
elseif strcmp(pref,'right')
    x0 = max(x0_vec);
else
    x0 = x0_vec;
end

end