function [x,t,cf]=SimCalCf(y0)

smtX=0.005;

%indexStart=128;
%indexStart=20;
indexStart=140;

%------- Find peak at each grid point.
% [~,indexT]=max(y0.Sxy(indexStart:end,:),[],1);
% t=y0.t(indexStart:end);
% x=y0.x;
% t=t(indexT);

% %----Exclude dt=0
% 
% dtTmp=diff(t);
% t=t(dtTmp~=0);
% dt=diff(t);
% x=x(dtTmp~=0);
% dx=diff(x);
% 
%   cf=dx./dt;
%   x=x(2:end);


% %------Find peak at each time step

[~,indexX]=max(y0.Sxy(indexStart:end,:),[],2);
t=y0.t(indexStart:end);
x=y0.x(indexX);
% %--Exclude dx=0
dxTmp=diff(x);
x=x(dxTmp~=0);
t=t(dxTmp~=0);

%-------------David's averaging method
for j=1:length(x)
    index1=find(x-x(j)<-smtX/2,1,'last');
    if(isempty(index1))
        index1=1;
    end
    index2=find(x-x(j)>smtX/2,1,'first');
    
    if(isempty(index2))
        index2=length(x);
    end
    cf(j)=(x(index1)-x(index2))./(t(index1)-t(index2));
end


