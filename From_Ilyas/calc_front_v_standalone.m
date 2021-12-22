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

end