function [s1 s2 angle]=calcPrincipalStresses(acq)


a=(acq.Sxx+acq.Syy)/2;
b=(((acq.Sxx-acq.Syy)/2).^2+acq.Sxy.^2).^0.5;

s1=a+b;
s2=a-b;
angle=-0.5*atan(2*acq.Sxy./(acq.Sxx-acq.Syy));%compression is positive
angle=angle*180/pi;