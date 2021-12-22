function T=FTIR(phi1,d)
%phi1 angle of incident deg
%D [\lambda]
%Frustrated total internal reflection: A demonstration and review S. Zhu
%et.al 1985

phi1=phi1/180*pi;

n1=1.5;
n2=1;
n3=1.5;
n=n3/n1;
N=n1/n2;

lambda=1;
l=lambda/(4*pi)./(n1^2*sin(phi1).^2-n2^2).^0.5;
y=d./(2*l);
%y=2*pi*d/lambda*(n1^2*sin(phi1).^2-n2^2).^0.5;

alphaPerp=(N^2-1)*(n^2*N^2-1)./(4*N^2*cos(phi1).*(N^2*sin(phi1).^2-1).*(n^2-sin(phi1).^2).^0.5);
betaPerp=((n^2-sin(phi1).^2).^0.5+cos(phi1)).^2./(4*cos(phi1).*(n^2-sin(phi1).^2).^0.5);

T=1./(alphaPerp.*sinh(y).^2+betaPerp);