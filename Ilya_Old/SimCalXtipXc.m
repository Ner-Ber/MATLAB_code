function [anl]=SimCalXtipXc(y0)



%------Find peak at each time step
indexT1=400;
indexT2=length(y0.Sxy(:,1))-10;

Sxy=y0.Sxy(indexT1:indexT2,:);
tTip=y0.t(indexT1:indexT2);

[~,indexX]=max(Sxy,[],2);
x=y0.x;
xTip=y0.x(indexX);

tau_r=3.704e6;
%tau_r=0;

for j=1:length(tTip)
    SxyTmp=Sxy(j,:);
    indexR=find(SxyTmp>tau_r,1);%find end of the cohesive zone
    xCend(j)=x(indexR-1)+( x(indexR)-x(indexR-1) )*(1-( SxyTmp(indexR)-SxyTmp(indexR-1))/0.5e6);%some interpolation for better resolution
end

anl.tTip=tTip;
anl.xTip=xTip;
anl.xCend=xCend';
anl.Cf=diff(xTip)./diff(tTip);