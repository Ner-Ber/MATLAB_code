function outmu=ds_v2mu_1990(involts)
% data is mm vs. volts
ds1990ns=dlmread('d:\frics\dscal\ds1990_nearside.txt');
%the 1000 is to convert to microns
outmu=1000*interp1(ds1990ns(:,2),ds1990ns(:,1),involts,'pchip')';