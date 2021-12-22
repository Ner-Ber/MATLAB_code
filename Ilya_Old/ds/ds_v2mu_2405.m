function outmu=ds_v2mu_2405(involts)
% data is mm vs. volts
ds2405fs=dlmread('c:\frics\dscal\ds2405_farside.txt');
%the 1000 is to convert to microns
outmu=1000*interp1(ds2405fs(:,2),ds2405fs(:,1),involts,'pchip');