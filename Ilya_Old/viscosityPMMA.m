function[t,Uxy]=viscosityPMMA
% the function converts a profile of stress to strain according to the
% dynamic properties given by dynamicModulusPMMA

%--------------define stress profile
t=(-7e-2:1e-5:7e-2);
Sxy=(1-tanh(t/30e-6))*0.5*1e6;
[Sxy_w,f]=my_fft(Sxy,t(2)-t(1));

%--------convert to strain
mu_f=dynamicModulusPMMA(f);
Uxy_w=Sxy_w/2./mu_f;
[Uxy,~]=my_ifft(Uxy_w,f(2)-f(1));

%--plot results
subplot(1,3,1)
semilogx(f(f>0),real(mu_f(f>0)),'.-');
xlabel('f(Hz)')
ylabel('R{mu}')
hold all;
subplot(1,3,2)
semilogx(f(f>0),imag(mu_f(f>0))./real(mu_f(f>0)),'.-');
ylabel('I{mu}/R{mu}')
xlabel('f(Hz)')
hold all;
subplot(1,3,3)
plot(t*1e3,Uxy*1e3,'.-');
hold all;
xlabel('t(ms)');
ylabel('Uxy');
xlim([-10 50]);
%plot(t*1e3,(Sxy(1)/2/mu_inf-(Sxy(1)-Sxy)/2/mu_0)*1e3,'black.-');

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



