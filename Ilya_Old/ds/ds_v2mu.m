function outmu=ds_v2mu(involts,serialNum)
% data is mm vs. volts

path=['c:\frics\DScal\ds' num2str(serialNum) '_farside.txt'];

if exist(path)
        
dsFs=dlmread(path);
%if signal's max signal is larger then the calibration -> renormalize.
dsVMax=max(dsFs(:,2));
involtsVMax=max(involts);
if(dsVMax<involtsVMax)
    involts=involts/involtsVMax*dsVMax;
end

%the 1000 is to convert to microns
outmu=interp1(dsFs(:,2),dsFs(:,1),involts,'pchip')*1000;

else
    outmu=involts;
end