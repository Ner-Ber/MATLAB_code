%-- load
[array,wavelength,aspect,pxlsize] = ReadOPD('C:\Users\NB\Google Drive\JAY_lab\Upper block scratches\elsa block with pattern\Sample_4.OPD',nan,0);
DX = 371:479;
DY = 1:480;
%-- crop
portion = array(DY,DX);
meanPort = mean(portion,2);
%-- fit linear
XX = (1:length(meanPort))';
meanPort = meanPort(:);
PP = polyfit(XX(~isnan(meanPort)),meanPort(~isnan(meanPort)),1);
YY = polyval(PP,XX);
arrayFixed = bsxfun(@minus,array,YY(:))./cos(atan(PP(1)));