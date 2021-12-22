% simulates propagating back from the camera 
lambda = 633*10^-9;
% z_prop =-2;
dz = 0.01; % m
dx = 6.4815e-06; % 6.5 micro meter pixel
x = 0:dx:dx*1400;
x = x - 700*dx;
E_in = reconstE1(201:1601,:); % to make x and y equal in size

figure;
z = zeros(1401,1401,40);
for j=1:150
    j
z_prop = j*dz;
E_out=Ang_spec_prop_2d(E_in,x,z_prop,lambda);
subplot(1,2,1);imagesc(angle(E_out));
title(strcat('phase.   z = -',num2str(z_prop),'m'));
subplot(1,2,2);imagesc(abs(E_out).^2);
title(strcat('amplitude^2.   z = -',num2str(z_prop),'m'));
pause(0.005)
% z(:,:,j) = abs(E_out).^2;
end
% 
% [X Y]=meshgrid(x);
% Z=zeros(1401);
% for i = 1:1401
%     for j = 1:1401
%     Z(i,j) = z(i,j,40);
%     end
% end
% 
% figure;surf(X,Y,Z,'EdgeColor','none');