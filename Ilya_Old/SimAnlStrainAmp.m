function anl=SimAnlStrainAmp(y0,y)

xStart=20e-3;
xEnd=220e-3;

[xf,cf]=SimCalCf(y0);
cf=smooth(cf,31);

[~,indexTStart]=min(abs(xf-xStart));
[~,indexTEnd]=min(abs(xf-xEnd));
[~,indexXStart]=min(abs(y.x-xStart));
[~,indexXEnd]=min(abs(y.x-xEnd));

anl.cf=cf(indexTStart:indexTEnd);

% %------on the interface
% Uxx0=0.5*diff(y0.slip,1,2)./repmat(diff(y0.x)',length(y0.slip(:,1)),1);
% Uxx0=Uxx0(indexTStart:indexTEnd,indexXStart:indexXEnd);
% anl.Uxx0M=-min(Uxx0'-Uxx0(1,end));


%---------off fault
Uyy=y.Uyy(indexTStart:indexTEnd,indexXStart:indexXEnd);
Uxx=y.Uxx(indexTStart:indexTEnd,indexXStart:indexXEnd);
anl.UyyM=max(Uyy'-Uyy(1,end));
anl.UxxM=-min(Uxx'-Uxx(1,end));

figure(20);
subplot(1,2,1);hold all;
plot(anl.cf,anl.UxxM*1e3,'.-');
ylabel('Uxx(mStrain)') 
subplot(1,2,2);hold all;
plot(anl.cf,anl.UyyM*1e3,'.-');
ylabel('Uyy(mStrain)') 