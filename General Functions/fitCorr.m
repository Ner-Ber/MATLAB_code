function [f,xx,yy]=fitCorr(r,lags,pointsNum)

m=find(r==max(r));
r2fit=r(m-pointsNum:m+pointsNum);
lags2fit=lags(m-pointsNum:m+pointsNum);
f=fit(lags2fit,r2fit,'poly2');
% int=integrate(f,lags2fit,0);

xx=lags2fit(1):0.0001:lags2fit(end);
yy=f.p1*xx.^2 + f.p2*xx + f.p3;

% figure;hold all
% plot(f,lags2fit,r2fit,'x')
% plot(xx,yy)
end