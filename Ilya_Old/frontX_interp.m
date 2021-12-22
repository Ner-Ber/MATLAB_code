function x=frontX_interp(frontX,frontT,p,t)


%x=csaps(frontT,frontX,p,t);
s=csaps(frontT,frontX,p);
s = fnxtr(s);
x=fnval(s,t);


% if (acqE.intVdt(end)<0)
%  [vdtMax,index]=findpeaks(acqE.intVdt);%Sometimes interpolation is not monotonic
%  acqE.intVdt(index(end):end)=vdtMax(end);
% end
% 
% if acqE.intVdt(1)>200
%  [vdtMin,index]=findpeaks(-acqE.intVdt);%Some interpolation is not monotonic
%  acqE.intVdt(1:index(1))=-vdtMin(1);
% end
% 
