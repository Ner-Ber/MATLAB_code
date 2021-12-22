function [UxxO,UxyO,UyyO]=calc_shear_sensitivity3(Uxx,Uxy,Uyy,sg_angle)
%Corection for transvers and shear sensitivity for general sg_angle.
%The input strain should not be rotated.
g=0.15; %shear sensititvity
g2=0;
k=0;

% %----- sg_angle=0
% %------Uxy
UxyO=1/(1-k)*Uxy;
%UxxO=1/(1+g-k*(k-g))*(Uxx+(g-k)*Uyy-g*(1+k+g)/(1-k)*Uxy );
UxxO=1/(1+g-k*(k-g))*(Uxx+(g-k)*Uyy-g2*(1+k+g)/(1-k)*Uxy ); 
UyyO=Uyy-k*UxxO+g2/(1-k)*Uxy;

%------Rotate strain
[UxxO,UxyO,UyyO]=RotateStrain(UxxO,UxyO,UyyO,sg_angle);


function [Uxx,Uxy,Uyy]=RotateStrain(Uxx,Uxy,Uyy,sg_angle)

%------------------- rotate the strain (correction)
for j=1:length(Uxx(1,:))
    u(1,1,:)=Uxx(:,j);
    u(1,2,:)=Uxy(:,j);
    u(2,1,:)=Uxy(:,j);
    u(2,2,:)=Uyy(:,j);
    R=[cos(sg_angle(j)) sin(sg_angle(j)); -sin(sg_angle(j)) cos(sg_angle(j))]; % rotation => u=R(theta)u'R(-theta)
    u=multiprod(u,R);
    R=[cos(sg_angle(j)) -sin(sg_angle(j)); sin(sg_angle(j)) cos(sg_angle(j))];
    u=multiprod(R,u);
    Uxx(:,j)=u(1,1,:);
    Uxy(:,j)=u(1,2,:);
    Uyy(:,j)=u(2,2,:);
end
