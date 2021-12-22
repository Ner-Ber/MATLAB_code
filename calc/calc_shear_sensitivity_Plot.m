function calc_shear_sensitivity_Plot(acqS,gV,color)
%acqS should be not rotated
%gV=[0,0.1,0.95,-0.08];

% gV(1)=-0.03;
% gV(2)=0.92;
% gV(3)=-0.22;
%[U1,U2,U3]=calcSgStrain(acqS.Uxx,acqS.Uxy,acqS.Uyy);
U1=acqS.U1;
U2=acqS.U2;
U3=acqS.U3;

[U1O,U2O,U3O]=calc_shear_sensitivity4(U1,U2,U3,gV);
[~,~,~,Uxx,Uyy,Uxy,~]=calculate_stress_strain(U1O,U2O,U3O,acqS.sg_angle);

figure;
hold all;

%index=[6,7,8,14,15,16];
index=1:4;

for j=1:length(index)
plot(U1O(:,index(j)),'-');
my_legend_add([num2str(acqS.x_sg(index(j))),' , U1']);legend off;
plot(U2O(:,index(j)),'-');
my_legend_add('U2');legend off;
plot(U3O(:,index(j)),'-');
my_legend_add('U3');legend off;
plot(U1O(:,index(j))+U3O(:,index(j))-U2O(:,index(j)),'-');
my_legend_add('U1+U3-U2');legend off;
end
title(['gV=' num2str(gV)]);
% a=get(gca,'Children');
% for j=1:length(a) set(a(end-j+1),'color',color{j}); end
%xlim([0 12e4]);
