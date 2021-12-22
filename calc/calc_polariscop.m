function I=calc_polariscop(acq,polAngle)
%I=calc_photoElasticity(acq,polAngle)
%polAngle is in degree
C=60E-13;
L=6E-3;
Lambda=670E-9;

[s1 s2 angle]=calcPrincipalStresses(acq);
angle=(angle-polAngle)*pi/180;
R=2*pi/Lambda*L*C*(s1-s2)*1E6;


I.I1=(sin(2*angle).^2);
I.I2=sin(R/2).^2;
I.I=I.I1.*I.I2;
I.t=acq.t;
