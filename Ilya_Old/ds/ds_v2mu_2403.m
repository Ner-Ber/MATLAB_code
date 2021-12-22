function outmu=ds_v2mu_2403(involts)
% data is mm vs. volts
ds2403fs=dlmread('c:\frics\dscal\ds2403_farside.txt');
%the 1000 is to convert to microns
outmu=interp1(ds2403fs(:,2),ds2403fs(:,1),involts,'pchip')*1000;