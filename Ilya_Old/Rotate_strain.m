function [Uxx,Uxy,Uyy]=Rotate_strain(Uxx,Uxy,Uyy,sg_angle)

if (length(sg_angle)==1)
sg_angle=sg_angle*ones(1,length(Uxx(1,:)));
end

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