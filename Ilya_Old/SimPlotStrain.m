function SimPlotStrain(y,indexT,figNum)

x0=0;

figure(figNum(1)); 
plot((y.x-x0)*1e3,(y.Uxy(indexT,:))*1e3-0.87,'.-');
hold all;
figure(figNum(2)); 
plot((y.x-x0)*1e3,(y.Uxx(indexT,:)-y.Uxx(1,end))*1e3,'.-');
hold all;
figure(figNum(3)); 
plot((y.x-x0)*1e3,(y.Uyy(indexT,:)-y.Uyy(1,end))*1e3,'.-');
hold all;


%on the interface
%[~,index]=max(y0.Sxy');
%[~,indexT]=min(abs(y0.x(index)-0.08));
%----plot Sxy
%plot((y0.x-y0.x(index(indexT)))*1e3,y0.Sxy(indexT,:)*1e-6-0.74*5,'.-')
%-----plot Uxx
%plot((y0.x(2:end)-y0.x(index(indexT)))*1e3,0.5*diff(y0.slip(indexT,:))'./diff(y0.x)*1e3,'.-');