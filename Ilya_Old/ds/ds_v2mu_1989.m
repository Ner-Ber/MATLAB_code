function outmu=ds_v2mu_1989(involts)
% data is mm vs. volts
ds1989fs=dlmread('d:\frics\dscal\ds1989_farside.txt');
%the 1000 is to convert to microns
outmu=interp1(ds1989fs(:,2),ds1989fs(:,1),involts,'pchip');