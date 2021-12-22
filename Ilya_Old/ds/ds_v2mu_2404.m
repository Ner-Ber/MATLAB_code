function outmu=ds_v2mu_2404(involts)
% data is mm vs. volts
ds2404fs=dlmread('c:\frics\DScal\ds2404_farside.txt');
%the 1000 is to convert to microns
outmu=interp1(ds2404fs(:,2),ds2404fs(:,1),involts,'pchip')*1000;