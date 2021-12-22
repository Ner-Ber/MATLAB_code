function [Sxy_f,t_f,Sxy_s,t_s]=calc_dynamic_slow_stressFFT(Uxy_f,t_f,Uxy_s,t_s)
%The function get syncronized slow and fast measurements are calculates Sxy
%by fft

%--Combine slow and fast 
t=[t_s t_f];
Uxy=[Uxy_s Uxy_f];
[t,index]=sort(t);
Uxy=Uxy(index);

%-----interp to 1mus resolution
t_interp=(t(1):1e-6:t(end));
Uxy_interp=interp1(t,Uxy,t_interp);

Sxy_interp=calc_dynamic_stressFFT(Uxy_interp,1e-6);

Sxy_s=interp1(t_interp,Sxy_interp,t_s);
Sxy_f=interp1(t_interp,Sxy_interp,t_f);



