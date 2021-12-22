function [SqrtPref Ux] = Crack_UxMode2ThetaToZero(v,r,varargin)
% Ux = Crack_UxMode2ThetaToZero(v,r,theta)
% 
% Crack_UxMode2ThetaToZero will calculate the particle displacemnt Ux on
% the interface: theta=0. 
%
% INPUTS:
% [v]=[Cr]      front velocity in Cr units
% [r]=[m]       distance along interface in meters
% [theta]=[rad] angle from crack tip should be as small as possible (OPTIONAL)

theta = setDefaults4function(varargin,1e-12);

[Cd, Cs, Cr ,~, ~, ~, mu , G ,~]=CrackSolutionMaterialProperties;

k=Cs/Cd;%Broberg p.330
v=v*Cr;

alpha_d=(1-(v/Cd).^2).^0.5;
alpha_s=(1-(v/Cs).^2).^0.5;
D=4*alpha_d.*alpha_s-(1+alpha_s.^2).^2;

%-----Following Freund - p.234
% A=alpha_s.*v.^2./D/(1-nu)/Cs^2;% I think this is only for plain strain. For plane stress should be  A=(1+nu)*alpha_s.*v.^2./D/Cs^2;
%
% if(PlaneStrain==true)
%     K=(G*E/(1-nu^2)./A).^0.5; %plane strain 
% else
%     K=(G*E./A).^0.5; %plane stress
% end

%----- Following Broberg p.334,336
A=2*(1-k^2)*alpha_s.*v.^2./D/Cs^2;
K=(G*mu*4*(1-k^2)./A).^0.5;

if (theta>0)
    theta_d=atan(alpha_d.*tan(theta));
    index=find(theta_d<0);%matlab gives -pi/2<theta<pi/2, I want 0<theta<pi
    theta_d(index)=theta_d(index)+pi;
    
    theta_s=atan(alpha_s.*tan(theta));
    index=find(theta_s<0);%matlab gives -pi/2<theta<pi/2, I want 0<theta<pi
    theta_s(index)=theta_s(index)+pi;
else
    theta_d=atan(alpha_d.*tan(theta));
    index=find(theta_d>0);%matlab gives -pi/2<theta<pi/2, I want -pi<theta<0
    theta_d(index)=theta_d(index)-pi;
    
    theta_s=atan(alpha_s.*tan(theta));
    index=find(theta_s>0);%matlab gives -pi/2<theta<pi/2, I want -pi<theta<0
    theta_s(index)=theta_s(index)-pi;
    
end


Ux = -((2*alpha_s*K)./(mu*D*sqrt(2*pi)))*sqrt(r)*(2*sin(theta_d/2)-(1+alpha_s.^2)*sin(theta_s./2));



end