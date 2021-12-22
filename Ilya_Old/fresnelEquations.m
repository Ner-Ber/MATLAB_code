function [I E]=fresnelEquations(ni,nt,ai)
%ai in radians
if (nargin<3)
ai=0:pi/1000:pi/2;
end
at=asin(ni/nt*sin(ai));
c=nt*cos(at)./(ni*cos(ai));
E.ai=ai*180/pi;
I.ai=E.ai;

E.rPerp=-sin(ai-at)./sin(ai+at);%The notations are in reference to Plane of Incidence (not in reference to interface) s-polarization,TE
E.rParl=tan(ai-at)./tan(ai+at);% p-polarization,TM
E.tPerp=2*sin(at).*cos(ai)./sin(ai+at);
E.tParl=2*sin(at).*cos(ai)./sin(ai+at)./cos(ai-at);

I.TPerp=c.*abs(E.tPerp).^2;
I.TParl=c.*abs(E.tParl).^2;
I.RPerp=abs(E.rPerp).^2;
I.RParl=abs(E.rParl).^2;


figure;
subplot(2,2,[1 2]);
plot(I.ai,[I.TPerp ; I.TParl ; I.RPerp ; I.RParl],'.-');
title('intensity')
legend({'TPerp'  'TParl'  'RPerp' 'RParl'});legend('off')
subplot(2,2,3);
plot(E.ai,[abs(E.tPerp).^2;abs(E.tParl).^2;abs(E.rPerp).^2;abs(E.rParl).^2],'.-');
title('abs(E)^2');
subplot(2,2,4);
plot(E.ai,[angle(E.tPerp)-angle(E.tParl);angle(E.rPerp)-angle(E.rParl)]*180/pi,'.-');
title('phase difference');