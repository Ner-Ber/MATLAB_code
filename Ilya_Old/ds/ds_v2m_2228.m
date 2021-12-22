function outm=ds_v2m_2228(involts)
% data is mm vs. volts
ds2228fs=dlmread('d:\frics\dscal\ds2228_farside.txt');
%the 1/1000 is to convert to meters
outm=1/1000*interp1(ds2228fs(:,2),ds2228fs(:,1),involts,'pchip')';