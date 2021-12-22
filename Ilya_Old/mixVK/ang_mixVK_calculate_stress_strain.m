function [Sxx,Syy,Sxy,Uxx,Uyy,Uxy,note]=ang_mixVK_calculate_stress_strain(sg,exp_details)

sg_angle=exp_details.sg_angle;
angle=exp_details.angle;
sort_index=exp_details.axis_sort_index;

angle(1,:)=pi-angle(1,:);
angle(2,:)=pi/2+angle(2,:);

cosa=cos(angle);
sina=sin(angle);
%--------------------calculate strain stress
%-in the Camera axis frame (positive is oposite to pushing direction)
%%% a remettre
a1=repmat(cosa(1,:),length(sg(:,1)),1);
a2=repmat(cosa(2,:),length(sg(:,1)),1);
a3=repmat(cosa(3,:),length(sg(:,1)),1);
b1=repmat(sina(1,:),length(sg(:,1)),1);
b2=repmat(sina(2,:),length(sg(:,1)),1);
b3=repmat(sina(3,:),length(sg(:,1)),1);


f=( a2.*a3.^2.*b1.^2.*b2- ...
    a1.*a3.^2.*b1.*b2.^2 - ...
    a2.^2.*a3.*b1.^2.*b3+...
    a1.^2.*a3.*b2.^2.*b3+ ...
    a1.*a2.^2.*b1.*b3.^2 - ...
    a1.^2.*a2.*b2.*b3.^2 ).^(-1);

Uxy=0.5*f.*((-a3.^2.*b2.^2 + a2.^2.*b3.^2).*sg(:,1:3:end)+...
    (b3.^2.*b1.^2 - a1.^2.*b3.^2).*sg(:,2:3:end)+...
    (-a2.^2.*b1.^2 +a1.^2.*b2.^2).*sg(:,3:3:end));
Uxx=f.*((a3.*b2.^2.* b3 - a2.* b2.*b3.^2).* sg(:,1:3:end)+...
     (-a3.*b1.^2 .*b3 + a1 .*b1 .*b3.^2).* sg(:,2:3:end)+...
     (a2 .*b1.^2.* b2 - a1.* b1.* b2.^2).* sg(:,3:3:end));
Uyy=f.*((a2.*a3.^2.*b2 - a2.^2.*a3.*b3).*sg(:,1:3:end)+...
      (-a1.*a3.^2.*b1 + a1.^2.*a3.*b3).*sg(:,2:3:end)+...
      (a1.*a2.^2.*b1 - a1.^2.*a2.*b2).*sg(:,3:3:end));
%%% a remettre

% %%% temporary
% Uxy=sg(:,1:3:end);
% Uyy=sg(:,2:3:end);
% Uxx=sg(:,3:3:end);

% %temporary
% Uxy=(sg(:,1:3:end)-sg(:,3:3:end))/2;
% Uyy=sg(:,2:3:end);                     %calculated normal strain
% Uxx=(sg(:,1:3:end)+sg(:,3:3:end))-Uyy; %calculated transvers strain

Uxy(:,11:19)=-Uxy(:,11:19); %Vishay

% Uxx(:,1:10)=1.7*Uxx(:,1:10);
%------------sort the strains with ascending order of x_sg. x=0 is where the stoper is.
Uxy=Uxy(:,sort_index)*1000; %[mStrain]
Uyy=Uyy(:,sort_index)*1000; %[mStrain]
Uxx=Uxx(:,sort_index)*1000; %[mStrain]

%------------------- rotate the strain (correction)
% for j=1:length(Uxx(1,:))
%     u(1,1,:)=Uxx(:,j);
%     u(1,2,:)=Uxy(:,j);
%     u(2,1,:)=Uxy(:,j);
%     u(2,2,:)=Uyy(:,j);
%     R=[cos(sg_angle(j)) sin(sg_angle(j)); -sin(sg_angle(j)) cos(sg_angle(j))]; % rotation => u=R(theta)u'R(-theta)
%     u=multiprod(u,R);
%     R=[cos(sg_angle(j)) -sin(sg_angle(j)); sin(sg_angle(j)) cos(sg_angle(j))];
%     u=multiprod(R,u);
%     Uxx(:,j)=u(1,1,:);
%     Uxy(:,j)=u(1,2,:);
%     Uyy(:,j)=u(2,2,:);
% end

%---------now calc stresses using hook law and assuming Szz=0
% Cd=2700;
% Cs=1345;
% sigma=(2*(Cs/Cd).^2-1)/2/((Cs/Cd).^2-1);%plane stress

sigma=1/3;%poisson ratio
E=3E+9/10^6 ;%youngs modulus [MPa]
E=E/1000; %because U is in mStrain

%---------Plane Stress (Szz=0)
Sxx=-E/(1-sigma^2)*(Uxx+sigma*Uyy);  % "-" is added positive stress is compresion
Syy=-E/(1-sigma^2)*(Uyy+sigma*Uxx); % "-" is added positive stress is compresion
Sxy=E/(1+sigma)*Uxy;
note='Plane Stress';

%---------Plane Strain (Uzz=0)

% Sxx=-E/(1+sigma)/(1-2*sigma)*((1-sigma)*Uxx+sigma*Uyy); % "-" is added positive stress is compresion
% Syy=-E/(1+sigma)/(1-2*sigma)*((1-sigma)*Uyy+sigma*Uxx); % "-" is added positive stress is compresion
% Sxy=E/(1+sigma)*Uxy;
% note='Plane Strain';
% %----
% Uzz=-sigma/(1-sigma)*(Uxx+Uyy);
% Sxx=E/((1+sigma)*(1-2*sigma))*((1-sigma)*Uxx+sigma*(Uyy+Uzz));
% Syy=E/((1+sigma)*(1-2*sigma))*((1-sigma)*Uyy+sigma*(Uxx+Uzz));
% Sxy=E/(1+sigma)*Uxy;
