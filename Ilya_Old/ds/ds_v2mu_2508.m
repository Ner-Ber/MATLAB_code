function outmu=ds_v2mu_2508(involts)
% data is mm vs. volts
ds2508fs=dlmread('c:\frics\DScal\ds2508_farside.txt');
%if signal's max signal is larger then the calibration -> renormalize.
dsVMax=max(ds2508fs(:,2));
involtsVMax=max(involts);
if(dsVMax<involtsVMax)
    involts=involts/involtsVMax*dsVMax;
end

%the 1000 is to convert to microns
outmu=interp1(ds2508fs(:,2),ds2508fs(:,1),involts,'pchip')*1000;