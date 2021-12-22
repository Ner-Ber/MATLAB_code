function fraun_main()
a=1000 * 10^-9;       % distance from source to image plane
width = 10000 * 10^-9;
height= 10000 * 10^-9;
size = 700;
A0 = 3; %# amplitude of reference wave
A1 = 1; %# amplitude of object wave  A0>>A1: A0/A1>=3
lambda = 632 * 10^-9; % wavelength in nanometers

x=linspace(0,width,size); % vector from 0 to width
y=linspace(0,height,size); % vector from 0 to height
[X,Y]=meshgrid(x,y); % matrices of x and y values at each position

s=Super(A0, A1, X-(width/2), Y-(height/2), a, lambda); % size-by-size (700x700) 
r=rand(size); % 700x700 matrix of random values on [0 1]

im = zeros(size);
im(r<(s/(A0+A1))) = 1; %# do this all at once instead of pixel-by-pixel

% display the image
figure
imshow(im,[])
title('test image')

end % end of function Interference


% Super is now vectorized, so you can give it a matrix of values for x and y
function S = Super(refamp,objamp,x,y,a,lambda)
    r1 = sqrt(a.*a+x.*x+y.*y); % dot notation: multiply element-wise
    S = refamp+(objamp*cos(2*pi*r1/(lambda)));
end % end of function Super