% instructions for use:
% ==================
% from pdetool, after running the solver
% export mesh [p e t]
% export solution [u]

% close all
clc
x=-2.5:0.05:2.5;
z=-1:0.01:2;
xscale=2.8e-4; %[cm]
%defines the unit of length in the simulation, which is also the width of the sample

% epsilon=13;
slice=-0.28;
d=19.5e-7; %cm

uxy=tri2grid(p,t,u,x,z);%[volts]
[Ex,Ez]=gradient(-uxy,x,z);
Ex=Ex*1e-3/xscale;%[kV/cm]
Ez=Ez*1e-3/xscale;%[kV/cm]

% figure;
% imagesc(x,z,Ez)
% colormap(cool)
% colorbar
% set(gca,'YDir','Normal')
% title('Ez (kV/cm)/1V')
% axis image
% hold on
% line([x(1) x(end)],[slice, slice],'LineStyle','--','Linewidth',2)
% axes('position',[0.5,0.5,0.3,0.25])
% plot(x,-Ez(abs(z-slice)<0.01,:)) % 1e6 to transform from [kV] to [mV]

scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)/8 scrsz(4)/8 scrsz(3)*3/4 scrsz(4)*3/4])

axes1=subplot(4,1,[1 2]);
pdeplot(p,e,t,'xydata',u,'xystyle','off','contour','on','levels',20,'Colormap',jet,'Colorbar','on');
set(gca,'FontSize',16)
% axis image
box on
ylim([z(1) z(end)])
xlim([x(1) x(end)])
% xlabel('x (a.u.)')
ylabel('z (a.u.)')
title(['dipole=',num2str(1e7*d),'nm, ','z=',num2str(slice)])

hold on
line([x(1) x(end)],[slice, slice],'LineStyle','--','Linewidth',2)
pos1=get(axes1,'Position');
% quiver(x,z,Ex,Ez)

axes2=subplot(4,1,3);
pos2=get(axes2,'Position');
pos2(3)=pos1(3);
set(axes2,'Position',pos2)
set(gca,'FontSize',16)
%%%the energy is -eEz*d=[e]*[kV/cm]*[cm] so E*d is in keV units
plot(x,-abs(Ez(abs(z-slice)<eps,:))*d*1e6) % 1e6 to transform from [kV] to [mV]
axis tight
% plot(x,uxy(abs(z-slice)<eps,:)) % 1e6 to transform from [kV] to [mV]

ylabel('\DeltaE=-E_z\cdoted (meV)')
% xlabel('x (a.u.)')

axes3=subplot(4,1,4);
pos3=get(axes3,'Position');
pos3(3)=pos1(3);
set(axes3,'Position',pos3)
set(gca,'FontSize',16)
plot(x,Ex(abs(z-slice)<eps,:)) % 1e6 to transform from [kV] to [mV]
axis tight
ylabel('E_x (kV/cm)')
xlabel('x (a.u., 1=sample width)')
