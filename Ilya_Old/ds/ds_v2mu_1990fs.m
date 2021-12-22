function outmu=ds_v2mu_1990fs(involts)
% data is mm vs. volts
ds1990fs=dlmread('c:\frics\dscal\ds1990_farside.txt');
%the 1000 is to convert to microns
outmu=1000*interp1(ds1990fs(:,2),ds1990fs(:,1),involts,'pchip')';