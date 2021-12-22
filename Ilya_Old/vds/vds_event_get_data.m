function vds=vds_event_get_data(dirname,eventNum,startl,endl,smt)
%The function returns event structure (not normilized).
%vds.x,vds.t,vds.lines

%oded used those parameters as function arguments
totalLen=150;
interval=1;


%check max left and right limits
dirname=[dirname '\vds'];
[leftLim rightLim]=vdsCalcLims(dirname,1);
cd ..
eval(sprintf('cd %s',dirname));
ims=vdsReadImsSlice(eventNum,startl,interval,endl,smt);
if (ims(1).Ind == -1)
    return
end
vdsMeta = readVdsMeta;
%already smoothed each line
%trun from struct array of lines into a matrix
lines=[ims.Line]';
imsLen=length(ims);
if smt>1
%     x=totalLen-((1:1:length(lines(1,:)))-leftLim)./(rightLim-leftLim)*totalLen; %Ilya's system is mirror
        x=((1:1:length(lines(1,:)))-leftLim)./(rightLim-leftLim)*totalLen; % my system
    t=(-(vdsMeta.EventNumBuffers*vdsMeta.BufferNumImages)+vdsMeta.PostIms+(ims(1).Ind:interval:ims(end).Ind))*vdsMeta.FrameT;
       %not really first line anymore
   %firstLine=repmat(lines(normaLine,:),imsLen,1);
   % lines=lines./firstLine-1;
else
%     x=totalLen-((1:1:length(lines(1,:)))-leftLim)./(rightLim-leftLim)*totalLen; %Ilya's system is mirror
    x=((1:1:length(lines(1,:)))-leftLim)./(rightLim-leftLim)*totalLen;
    t=(-(vdsMeta.EventNumBuffers*vdsMeta.BufferNumImages)+vdsMeta.PostIms+(ims(1).Ind:interval:ims(end).Ind))*vdsMeta.FrameT;
    %not really first line anymore
%     firstLine=repmat(smooth(lines(normaLine,:),10)',length(ims),1);
end

t=t*10^3; %[msec]
vds.t=t';
vds.x=x;
vds.lines=lines;

cd ..
cd ..