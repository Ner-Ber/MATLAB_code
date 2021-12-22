function [Sxy]=calc_dynamic_stressFFT(Uxy,dt)
%v should be column oriented
%Uxy [mStrain]
%[dt] =s;

Uxy=Uxy*1e-3;%[Strain]
[Uxy_f,f]=my_fft(Uxy,dt);

mu_f=dynamicModulusPMMA(f);
Sxy_f=Uxy_f*2.*repmat(mu_f,1,size(Uxy_f,2));
[Sxy,~]=my_ifft(Sxy_f,f(2)-f(1));
Sxy=Sxy*1e-6;%convert to MPa


function[y,f]=my_fft(v,dt)

if length(v(1,:))>1 && length(v(:,1))==1
    v=v.';
end

L=length(v(:,1));
Fs=1/dt;
%NFFT = 2^nextpow2(L); % Next power of 2 from length of y
NFFT=L;%L-1;
Y = fft(v,NFFT);
f = Fs/2*linspace(0,1,NFFT/2+1);

if (rem(NFFT,2)==0)
f=[-f(end:-1:1) f(2:end-1)];
else
f=[-f(end:-1:1) f(2:end)];
end

y=fftshift(Y);


function[y,t]=my_ifft(v,df)

if length(v(1,:))>1 && length(v(:,1))==1
    v=v.';
end
v=ifftshift(v);

L=length(v(:,1));
T=df;

NFFT=L;%L-1;
y = ifft(v,NFFT);
t = T*linspace(0,1,L);

