% function E_out=Ang_spec_prop_2d(E_in,x,z_prop,lambda);
% propagates BACKWARD
% arguments:
% E_in: input field
% x - x spatial coordinate vector with real physical length units (assuming y is the same)
% z_prop - propagation distance
% lambda - wavelength (in same length units as x vector)

function E_out=Ang_spec_prop_2d(E_in,x,z_prop,lambda)
s=length(x);
% if mod(s,2)~=0
%     display('x vector length not even!')
%     pause
% end
temp=[-s/2:s/2-1]/s; % [-0.5 0.5] unitless
% xx=temp*resolution*s; % x vector in mm
kx=temp*(1/(max(x)-min(x)))/(temp(2)-temp(1))*2*pi; % fourier coordinates in 2pi/mm
alpha=(kx/2/pi)*lambda;  % unitless normalized fourier coordinates
[Alpha Betta]=meshgrid(alpha);
mu = 2*pi/lambda * sqrt(Alpha.^2+Betta.^2-1);

angs_prop_z=exp(mu*z_prop); % angular spectrum of propogation distance z


% imagesc(abs(angs_prop_z));
% colorbar
% if size(angs_prop_z,1)~=size(E_in,1)
%     angs_prop_z=transpose(angs_prop_z);
% end
E_out=ifft2(ifftshift(fftshift(fft2(E_in)).*angs_prop_z));