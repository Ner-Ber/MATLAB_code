% Program to calculate 1 dimensional Fresnel and Fraunhoffer diffraction
% patterns for single rectangular slit
% clear all
% close all
disp('THIS PROGRAM CALCULATES FRESNEL DIFFRACTION PATTERNS- VALID FOR ALL DISTANCE SCALES.WHERE FRESNEL NUMBER << 1 CAN USE FRAUNHOFFER APPROX')
% Define diffraction properties
wavelength=input('What is the wavelength of the light used to illuminate the diffracting aperture (nm)? ');
width=input('What is the aperture width (mm)? ');
distance=input('What is the distance from aperture to screen (m)? ');
wl=10^-9;
wdth=10^-3;
% Calculation of Fresnel number
F=((width*wdth)^2)/((wavelength*wl)*distance);
disp(['Fresnel number = ',num2str(F)])
% Determination of Fresnel or Fraunhoffer diffraction
if F<=0.1
 G=input(' Compute Fresnel (1) or Fraunhoffer case (2)? ')
else
 G=1;
end
%% %Fraunhoffer case
if G==2

 x=[-10:0.005:10]*wdth;
 y=abs(x)<=(width*wdth)/2;
 scaling1=(wavelength*wl)/((x(2)-x(1))*(length(x)-1));
 scaling2=(x(2)-x(1));
 x1=scaling2*x;
 endx=[-((length(x)-1)/2):((length(x)-1)/2)]*(wavelength*wl)/((length(x)-1)*scaling2);

 yscal=sum((y.^2)*(x(2)-x(1)))/sum((abs(fft(y)).^2)*(endx(2)-endx(1)));

 k=(2*pi)/(wavelength*wl);
 z=distance;

 h2=((exp(i*k*z))/(i*(wavelength*wl)*z))*(exp(((i*k)/(2*z))*(endx.^2)));

 figure(1) % Transmittance vs x

 plot(x,y)
 title(['Rectangular slit of width',num2str(width),'mm'])
 ylabel('Transmittance')
 xlabel('x position (m)')

 figure(2) % Unscaled

 F1=fft(y);
 FFS1=h2.*fftshift(F1);
 plot(endx,abs(FFS1*sqrt(yscal)));

 title(['Cross section of Fraunhoffer diffraction pattern produced by rectangular aperture of width',num2str(width),'mm(unscaled)'])
 ylabel('Relative amplitude')
 xlabel('Angular distribution (unscaled)')

 figure(9) % Rescaled

 plot(endx,abs(((FFS1).^2)*yscal))
 title(['Cross section of Fraunhoffer diffraction pattern produced by rectangular aperture of width',num2str(width),'mm(rescaled)'])
 ylabel('Intensity')
 xlabel('Angular distribution (radians)')

 figure(10)

 posendx=z*tan(endx);
 plot(posendx,abs(((FFS1).^2)*yscal))
 title(['Cross section of Fraunhoffer diffraction pattern produced by rectangular aperture of width',num2str(width),'mm(rescaled)'])
 ylabel('Intensity')
 xlabel('Spatial distribution (m)')



end
%% % Fresnel Case
if G==1
 tic
 x=[-10:0.000005:10]*wdth;
 y=abs(x)<=(width*wdth)/2;

 figure(4) % Transmittance vs x

 plot(x,y)
 title(['Rectangular slit of width',num2str(width),'mm'])
 ylabel('Transmittance')
 xlabel('x position')

 k=(2*pi)/(wavelength*wl);
 z=distance;





 scaling1=(wavelength*wl)*z/((x(2)-x(1))*(length(x)-1));
 scaling2=(x(2)-x(1));
 x1=scaling2*x;
 endx=[-((length(x)-1)/2):((length(x)-1)/2)]*scaling1;


 h2=((exp(i*k*z))/(sqrt(i*(wavelength*wl)*z)))*(exp(((i*k)/(2*z))*(endx.^2)));

 h1=((exp(((i*k)/(2*z))*x.^2)));
%
 fres11=(ifftshift(ifft(fft(y).*(h1)))); % Convolution attempt by computing kernel in fourier space- doesn't work!
 fres1=ifftshift(ifft((fft(y)).*(fft(h2)))); % Convolution kernel in real space- does work but lack of resolution in real space%

 yscal=sum(((y.*abs(h1)).^2)*(x(2)-x(1)))/sum((abs(fft(y.*h1)).^2)*(endx(2)-endx(1)));


 fres2=fftshift(fft(y.*h1))*sqrt(yscal).*h2;

 









%%
% h1=exp((i*pi*z*(wavelength*wl))*(x1.^2));
% figure(4) % Convolution attempt by using fourier space kernel
%
% title('Method 1')
%
% subplot(2,5,1)
%
% plot(x,y)
% title(['Rectangular slit of width',num2str(width),'mm'])
%
%
% subplot(2,5,6)
%
% plot(x,angle(y))
% title(['Phase over rectangular slit of width',num2str(width),' mm'])
%
%
% subplot(2,5,2)
%
% plot(x1,fftshift(abs(fft(y))))
% title('Fourier transform of aperture')
%
%
% subplot(2,5,7)
%
% plot(x1,unwrap(angle(fft(y))))
% title(['Phase of fourier transform of aperture'])
%
% subplot(2,5,3)
%
% plot(x,abs(h1))
% title(['PLot of Fresnel propagator'])
%
% subplot(2,5,8)
%
% plot(x,unwrap(angle(h1)))
% title('Phase of Fresnel propagator')
%
% subplot(2,5,4)
%
% plot(x,fftshift(abs(fft(h1))))
% title('Fourier transform of Fresnel propagator ignore!!!!')
%
% subplot(2,5,9)
%
% plot(x,unwrap(angle(fft(h1))))
% title('Phase of Fresnel propagator ignore¬!!!!')
%
%
% subplot(2,5,5)
%
%
% title('Convolution of aperture field and propagator')
% plot(endx,ifftshift(abs(ifft((fft(y)).*(h1)))))
%
% subplot(2,5,10)
%
% title('Phase of Convolution')
% plot(endx,unwrap(angle(ifftshift(ifft((fft(y)).*(h1))))))

%
%
%%
% h3=((exp(i*k*z))/(i*(wavelength*wl)*z))*(exp(((i*pi)/((wavelength*wl)*z))*endx.^2));
% h3=(exp(((i*pi)/((wavelength*wl)*z))*endx.^2));
% figure(5) % Kernel in real space
%
%
% subplot(2,5,1)
%
% plot(x,y)
% title(['Rectangular slit of width',num2str(width),'mm'])
%
%
% subplot(2,5,6)
%
% plot(x,angle(y))
% title(['Phase over rectangular slit of width',num2str(width),' mm'])
%
%
% subplot(2,5,2)
%
% plot(x1,fftshift(abs(fft(y))))
% title('Fourier transform of aperture')
%
%
% subplot(2,5,7)
%
% plot(x1,unwrap(angle(fft(y))))
% title(['Phase of fourier transform of aperture'])
%
% subplot(2,5,3)
%
% plot(x,abs(h3))
% title(['PLot of Fresnel propagator'])
%
% subplot(2,5,8)
%
% plot(x,unwrap(angle(h3)))
% title('Phase of Fresnel propagator')
%
% subplot(2,5,4)
%
% plot(x,fftshift(abs(fft(h3))))
% title('Fourier transform of Fresnel propagator ignore!!!!')
%
% subplot(2,5,9)
%
% plot(x,angle(fft(h3)))
% title('Phase of Fresnel propagator ignore¬!!!!')
%
%
% subplot(2,5,5)
%
%
% title('Convolution of aperture field and propagator')
% plot(endx,ifftshift(abs(ifft((fft(y)).*(h3)))))
%
% subplot(2,5,10)
%
% title('Phase of Convolution')
% plot(endx,angle(ifftshift(ifft((fft(y)).*(h3)))))




%%
 




 figure(6) % Fourier transform of product of aperture and propagator

 title('Method 2')

 subplot(2,4,1)

 plot(x,y)
 title(['Rectangular slit of width',num2str(width),'mm'])

 subplot(2,4,5)

 plot(x,angle(y))
 title('Phase of rect function')

 subplot(2,4,2)

 plot(x,abs(h1))
 title('Fresnel propagator')

 subplot(2,4,6)

 plot(x,unwrap(angle(h1)))
 title('Phase of Fresnel propagator')

 subplot(2,4,3)

 plot(abs(y.*h1))
 title('Aperture x propagator')

 subplot(2,4,7)

 plot(unwrap(angle(y.*h1)))
 title('Phase of aperture x propagator')

 subplot(2,4,4)

 plot(abs(fftshift(fft(y.*h1))))
 title('Fourier transform of aperture x propagator')

 subplot(2,4,8)

 plot(unwrap(angle(fftshift(fft(y.*h1)))))
 title('Phase of fourier transform of aperture x propagator')

 figure(7) % Rescaled intensity plot


 plot(endx,abs(fres2).^2)
 ylabel('Intensity')
 xlabel('x position (m)')

 figure(8)

 angleendx=atan(endx/z);

 plot(angleendx,nonzeros(abs(fres2).^2))
 ylabel('Intensity')
 xlabel(['Angular distribution at distance ',num2str(z),'m (rad)'])




% subplot(2,5,9)
%
% plot(angle(h2))
%
% subplot(2,5,5)
% 
% 42
% plot((abs((fft(y).*h1))))
%
% subplot(2,5,10)
%
% plot(angle(fftshift(fft(y).*h1)))


% figure(5) % Amplitude spectrum
%
% subplot(1,2,1)
% plot(endx,absfres2)
% title(['Cross section of Fresnel diffraction pattern produced by rectangular aperture of width',num2str(width),'mm
% (rescaled) method 2'])
% ylabel('Relative amplitude')
% xlabel('x position (m)')
% disp(['Fresnel number = ',num2str(F)])
%
% subplot(1,2,2)
% plot(x,absfres1)
% title(['Cross section of Fresnel diffraction pattern produced by rectangular aperture of width',num2str(width),'mm
% (rescaled) method 1'])
% ylabel('Relative amplitude')
% xlabel('x position (m)')
%
% figure(6) % Intensity spectrum
%
% subplot(1,2,1)
% plot(endx,(absfres2).^2)
%
% title(['Cross section of Fresnel diffraction pattern produced by rectangular aperture of width',num2str(width),'mm
% (rescaled) method 2'])
% ylabel('Relative intensity')
% xlabel('x position (m)')
%
%
% subplot(1,2,2)
% plot(endx,(absfres1).^2)
%
% title(['Cross section of Fresnel diffraction pattern produced by rectangular aperture of width',num2str(width),'mm
% (rescaled) method 2'])
% ylabel('Relative intensity')
% xlabel('x position (m)')
%
 t=toc

end