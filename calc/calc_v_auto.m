function [x v]=calc_v_auto(exper,event,Xstart,Xend,tStart,tEnd)

%Also works for rupture propagating in negative direction. For bi directional propagation use seperatly.

exp_details=expDetailsRead(exper);
lines=2:7;
smt=15;
smtX=5;
%-------Find tStart and tEnd
if (nargin==4)
phE=phantomGetLines(exper,event,-6,'min',6,1000,3,'',lines);    
A=get_A_at_x(phE,[Xstart, Xend],5);
index=find(subtruct_norm(A.lines(:,1),1)<0.95,1,'first');
tStart=A.t(index);
index=find(subtruct_norm(A.lines(:,2),1)<0.95,1,'first');
tEnd=A.t(index);
end

phE=phantomGetLines(exper,event,tStart-1,'min',tEnd+2,1000,smt,'',lines);
phE.firstLine=mean(phE.lines(1:20,:),1);


[frontRaw front]=calc_frontXT_from_XAN(phE,tStart,tEnd,7,Xstart,Xend);

%-calc front v (set the smooth parameters at the function

[x v]=calc_front_v(frontRaw.x,frontRaw.t,smtX);
%[x v]=calc_front_v(front.x1,front.t);




function [xv,v]=calc_front_v(x,t,smtX)

xv=x(2:end)-0.5*(mean(diff(x)));

% if(x(1)<x(end))
%     [~,index]=min(abs(x-50));
%     [~,index2]=min(abs(x-90));
% else
%     [~,index]=min(abs(x-120));
%     [~,index2]=min(abs(x-80));
% end
%v=smooth(x(2:end),diff(x)./diff(t),5)';

v=diff(x)./diff(t);
v=my_smooth2(xv,v,smtX);



%vtmp=smooth(diff(x)./diff(t),3)';
% v1=vtmp(1:index);
% vtmp=smooth(diff(x)./diff(t),3)';
% v2=vtmp(index+1:index2);
% vtmp=smooth(diff(x)./diff(t),3)';
% v3=vtmp(index2+1:end);
% v=[v1 v2 v3];

