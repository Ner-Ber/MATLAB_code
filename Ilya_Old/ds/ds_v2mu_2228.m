function outmu=ds_v2mu_2228(involts)
% data is mm vs. volts
ds2228fs=dlmread('d:\frics\dscal\ds2228_farside.txt');
%the 1000 is to convert to microns
outmu=1000*interp1(ds2228fs(:,2),ds2228fs(:,1),involts,'pchip');