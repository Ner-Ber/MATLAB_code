function [e]=anls_ds(exper,event,Xstart,Xend)

xds=85;
%Also works for rupture propagating in negative direction. For bi directional propagation use seperatly.

%-------Find tStart and tEnd
if (nargin==4)
phE=phantomGetLines(exper,event,-2,'min',1,1000,3,'',2:7);    
A=get_A_at_x(phE,[Xstart, Xend],5);
index=find(subtruct_norm(A.lines(:,1),1)<0.9,1,'first');
tStart=A.t(index);
index=find(subtruct_norm(A.lines(:,2),1)<0.9,1,'first');
tEnd=A.t(index);
end

phE=phantomGetLines(exper,event,tStart-1,'min',tEnd+2,1000,3,'',3:6);
phE.firstLine=mean(phE.lines(1:20,:),1);
acqE=acq132_event_get_data(exper,event,'start','end',1,'ds2511');


%----detecte and smooth the front

[frontRaw front]=calc_frontXT_from_XAN(phE,tStart,tEnd,3,Xstart,Xend);

%----------
x=frontRaw.x; %another option is by front.x1
t=frontRaw.t;

[~,index]=min(abs(x-xds));
e.tds=t(index);
[~,index]=min(abs(acqE.t-e.tds));
e.ds0=mean(acqE.ds2511(index-100:index-50));
e.xds=xds;
e.acqE=acqE;




