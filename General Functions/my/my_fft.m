function[y,f]=my_fft(v,dt)

if length(v(1,:))>1 && length(v(:,1))==1
    v=v';
end

L=length(v(:,1));
Fs=1/dt;
%NFFT = 2^nextpow2(L); % Next power of 2 from length of y
NFFT=L-1;
Y = fft(v,NFFT);
f = Fs/2*linspace(0,1,NFFT/2+1);
y=abs(Y(1:NFFT/2+1,:)); %single sided ampiltude spectrum


% Plot single-sided amplitude spectrum.
%set(0,'DefaultFigureWindowStyle','docked') ;

% plot(f,y,'.-');
% figure;
% loglog(f,y,'.-') 
% title('Single-Sided Amplitude Spectrum of y(t)')
% xlabel('Frequency (Hz)')
% ylabel('|Y(f)|')
% 
