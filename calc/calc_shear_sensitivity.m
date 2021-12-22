function [UxxO,UxyO,UyyO]=calc_shear_sensitivity(Uxx,Uxy,Uyy,sg_angle)
%Corection for transvers and shear sensitivity for general sg_angle.
%The input strain should not be rotated.
g=0.15; %shear sensititvity
k=0;

% %----- sg_angle=0
% %------Uxy
UxyO=1/(1-k)*Uxy;

%--------Uxx
Uxx_Uyy_Sign=sign(  (1-k)*Uxx-(1+k)*Uyy+2*g*(1-k^2)/(1-k)*abs(Uxy)  ); %terms with k*g are neglected.

UxxOP=1/(1+g-k*(k-g))*( Uxx+(g-k)*Uyy+g*(1+k-g)/(1-k)*abs(Uxy) ); %--------For UxxO-UyyO>0
UxxOM=1/(1-g-k*(k+g))*( Uxx+(-g-k)*Uyy+g*(1+k+g)/(1-k)*abs(Uxy) ); %--------For UxxO-UyyO<0

UxxO=0*UxxOP;
UxxO(Uxx_Uyy_Sign>0)=UxxOP(Uxx_Uyy_Sign>0);
UxxO(Uxx_Uyy_Sign<0)=UxxOM(Uxx_Uyy_Sign<0);

%-------Uyy
UyyO=Uyy-k*UxxO-g/(1-k)*abs(Uxy);

 [U1,U2,U3]=calcSgStrain(Uxx,Uxy,Uyy);
%   UyyO=Uyy;
%   UxxO=Uxx;
%   UxyO=Uxy;

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
