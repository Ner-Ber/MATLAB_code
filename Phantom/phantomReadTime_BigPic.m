function [startl,interval,endl,t]=phantomReadTime_BigPic(exp_dir,eventNum,Tstart,Tinterval,Tend,timeBase)
%eventNum=0 is slow acquisition. timeBase=1000 -> msec
 %this function is responsible for all the time logic
 
if nargin<6
    timeBase=1;
end
%--------Get triggerDelay 
if ( exist([exp_dir '\exp_details.txt'],'file') && (eventNum~=0))%for slow no trigger delay
    expDetails=expDetailsRead(exp_dir);
    triggerDelay=expDetails.triggerDelay;
else
    triggerDelay=0;
end

%--------Time Logic

fid= fopen([exp_dir '\PhBig\' num2str(eventNum) 'Time.bin'],'r');
t=(fread(fid,inf,'double')+triggerDelay)*timeBase;
fclose(fid);

if(strcmp('end',Tend)||Tend>max(t))
    endl=length(t);
elseif(Tend<min(t))
    error('Tend out of limits');
else
    [~,endl]=min(abs(t-Tend));
end

if(strcmp('start',Tstart)||Tstart<min(t))
    startl=1;
elseif(Tstart>max(t))
    error('Tstart out of limits');
else
    [~,startl]=min(abs(t-Tstart));
end

if (strcmp(Tinterval,'min'))
    interval=1;
else
    interval=ceil(Tinterval/mean(diff(t(1:20))));
end

t=t(startl:interval:endl);
