function gVOut=calc_shear_sensitivity_fit(acq)
%acq should not be rotated

%sgIndex=[6,7,8,14,15,16];
%[U1,U2,U3]=calcSgStrain(acq.Uxx(5e3:10:3e5,sgIndex),acq.Uxy(5e3:10:3e5,sgIndex),acq.Uyy(5e3:10:3e5,sgIndex));
%sg_angle=acq.sg_angle(sgIndex);

%[U1,U2,U3]=calcSgStrain(acq.Uxx,acq.Uxy,acq.Uyy);
%sg_angle=acq.sg_angle;

U1=acq.U1;
U2=acq.U2;
U3=acq.U3;
sg_angle=acq.sg_angle;

%permut=[16 13:15];
%permut=[4,1,2,3];
permut=1:length(acq.U1);
 permutRand=randperm(length(permut));
 permut=permut(permutRand);
 

U1=U1(:,permut);
U2=U2(:,permut);
U3=U3(:,permut);
sg_angle=acq.sg_angle(permut);

gV0= [0,0]; %First check what gV0 stands for in lsq
lsqM=@(gV) lsq(U1,U2,U3,sg_angle,gV);
gVOut = lsqnonlin( lsqM ,gV0);

%%Plot for check
% figure;
% %---no correction
% [~,~,~,Uxx,Uyy,Uxy,~]=calculate_stress_strain(U1,U2,U3,sg_angle);
% %trace=Uxx+Uyy;
% %f=[diff(trace,1,2)];
% f=[diff(Uxx,1,2);diff(Uxy,1,2);diff(Uyy,1,2)];
% f=f(:);
% plot(f);

% %--after lsq
% hold all
% f=lsqM(gVOut);
% plot(f);
% title(num2str(gVOut));

function f=lsq(U1,U2,U3,sg_angle,gV)

%gV=[gV(1) 0.1 gV(2:3)];
%gV=[0 0.1 gV];
gV=[0 gV(1) 0.95 gV(2)];
%gV=[0 gV 1 0]; %for trace fit

[U1tmp,U2tmp,U3tmp]=calc_shear_sensitivity4(U1,U2,U3,gV);
[~,~,~,Uxx,Uyy,Uxy,~]=calculate_stress_strain(U1tmp,U2tmp,U3tmp,sg_angle);

% %-----Fit trace
%trace=Uxx+Uyy;
%f=[diff(trace,1,2)];
%f=trace--3.3e-3;

f=[diff(Uxx,1,2);diff(Uxy,1,2);diff(Uyy,1,2)];
%  f=[diff(Uxx,1,2);diff(Uxy,1,2);diff(Uyy,1,2);Uxx(1:end-1)+0.61*Uyy(1:end-1)];
% f=[diff(Uxx,1,2);diff(Uyy,1,2);Uxx(1:end-1)+0.61*Uyy(1:end-1)];

% 
%   UyyT=-3.3/3e6/pi/7.6e-3/50e-3*9.8/1.1;% [Strain/kg]
%    UxxT=2/3e6/pi/7.6e-3/50e-3*9.8/1.1;% [Strain/kg]
%    f=[Uxx-UxxT;Uyy-UyyT];

f=f(:);

