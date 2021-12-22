function [UxxO,UxyO,UyyO]=calc_shear_sensitivity2(Uxx,Uxy,Uyy,sg_angle)

%The resulting strain is calculated for sg_angle=0 and +-45 deg.

g=-0.15; %shear sensititvity
k=0;


% %----- sg_angle=0
% %------Uxy
%UxyO=1/(1-k)*Uxy;

% %--------Uxx
% Uxx_Uyy_Sign=sign(  (1-k)*Uxx-(1+k)*Uyy+2*g*(1-k^2)/(1-k)*abs(Uxy)  ); %terms with k*g are neglected.
% 
% UxxOP=1/(1+g-k*(k-g))*( Uxx+(g-k)*Uyy+g*(1+k-g)/(1-k)*abs(Uxy) ); %--------For UxxO-UyyO>0
% UxxOM=1/(1-g-k*(k+g))*( Uxx+(-g-k)*Uyy+g*(1+k+g)/(1-k)*abs(Uxy) ); %--------For UxxO-UyyO<0
% 
% UxxO=0*UxxOP;
% UxxO(Uxx_Uyy_Sign>0)=UxxOP(Uxx_Uyy_Sign>0);
% UxxO(Uxx_Uyy_Sign<0)=UxxOM(Uxx_Uyy_Sign<0);
% 
% %-------Uyy
% UyyO=Uyy-k*UxxO-g/(1-k)*abs(Uxy);
% 
% %------Rotate strain
% [UxxO,UxyO,UyyO]=RotateStrain(UxxO,UxyO,UyyO,sg_angle);


%---------sg_angle =-45 (with y axis)
% [U1,U2,U3]=calcSgStrain(Uxx,Uxy,Uyy);
% 
% Uxy=1/2*(U1+U3)-U2;
% Uxx=U3;
% Uyy=U1;
% 
% % %---for Uxy>0
%  UxyOP=1/(1-k+g)*(Uxy+g/2/(1-k)*abs(Uxx-Uyy));
% % %---for Uxy<0
%  UxyOM=1/(1-k-g)*(Uxy+g/2/(1-k)*abs(Uxx-Uyy));
%  UxyO=0*UxyOP;
%  UxyO(UxyOP>0)=UxyOP(UxyOP>0);
%  UxyO(UxyOP<0)=UxyOM(UxyOP<0);
%  
%  UxxO=1/(1-k^2)*(Uxx-k*Uyy-g*(1-k)*abs(UxyO));
%  UyyO=1/(1-k^2)*(Uyy-k*Uxx-g*(1-k)*abs(UxyO));

% %---------sg_angle =-45 (with y axis)

 [U1,U2,U3]=calcSgStrain(Uxx,Uxy,Uyy);
Uxy=-1/2*(U1+U3)+U2;
Uxx=U1;
Uyy=U3;

% %---for Uxy>0
 UxyOP=1/(1-k-g)*(Uxy-g/2/(1-k)*abs(Uxx-Uyy));
% %---for Uxy<0
 UxyOM=1/(1-k+g)*(Uxy-g/2/(1-k)*abs(Uxx-Uyy));

 % Assuming 1/(1-k+g) doesn't change the sign
 UxyO=0*UxyOP;
 UxyO(UxyOP>0)=UxyOP(UxyOP>0);
 UxyO(UxyOP<0)=UxyOM(UxyOP<0);
% 
 UxxO=1/(1-k^2)*(Uxx-k*Uyy-g*(1-k)*abs(UxyO));
 UyyO=1/(1-k^2)*(Uyy-k*Uxx-g*(1-k)*abs(UxyO));

function [Uxx,Uxy,Uyy]=RotateStrain(Uxx,Uxy,Uyy,sg_angle)

%------------------- rotate the strain
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
