function vdsPlotEventSliceXlim(dirname,eventNum,startl,endl,smt,baseFig)

%oded used those parameters as function arguments
totalLen=150;
interval=1;
normaLine=1;
timebase=1/1000; %[msec]

if nargin<6
    baseFig=gcf;
end
set(0,'DefaultFigureWindowStyle','docked')

dirname=[dirname '\vds'];

%return t vector for outside use

%check max left and right limits
[leftLim rightLim]=vdsCalcLims(dirname,1);
cd ..
eval(sprintf('cd %s',dirname));
ims=vdsReadImsSlice(eventNum,startl,interval,endl,smt);
if (ims(1).Ind == -1)
    return
end
vdsMeta = readVdsMeta;
cd ..
cd ..
% if ims(1).Im == -1
%     cd ..
%     return
% end
%already smoothed each line
%trun from struct array of lines into a matrix
lines=[ims.Line]';
imsLen=length(ims);
if smt>1
    %     x=totalLen-((1:1:length(lines(1,:)))-leftLim)./(rightLim-leftLim)*totalLen; %Ilya's system is mirror
    x=((1:1:length(lines(1,:)))-leftLim)./(rightLim-leftLim)*totalLen; % my system
    %         ; my system
    X=repmat(x,imsLen,1);
    t=(-(vdsMeta.EventNumBuffers*vdsMeta.BufferNumImages)+vdsMeta.PostIms+(ims(1).Ind:interval:ims(end).Ind))*vdsMeta.FrameT;
    T=repmat(t,length(x),1)'/timebase;
    %not really first line anymore
    firstLine=repmat(lines(normaLine,:),imsLen,1);
    %mesh(lines./firstLine);view(0,90);
    figure(baseFig);
    mesh(X,T,lines./firstLine);view(0,90);
    % mesh(lines./firstLine);view(0,90);
else
    %     x=totalLen-((1:1:length(lines(1,:)))-leftLim)./(rightLim-leftLim)*totalLen; %Ilya's system is mirror
    x=((1:1:length(lines(1,:)))-leftLim)./(rightLim-leftLim)*totalLen;
    X=repmat(x,length(ims),1);
    t=(-(vdsMeta.EventNumBuffers*vdsMeta.BufferNumImages)+vdsMeta.PostIms+(ims(1).Ind:interval:ims(end).Ind))*vdsMeta.FrameT;
    T=repmat(t,length(x),1)'/timebase;
    %not really first line anymore
         firstLine=repmat(smooth(lines(normaLine,:),10)',length(ims),1);
    %     mesh(lines);view(0,90);
    figure(baseFig);
    mesh(X,T,lines./firstLine);view(0,90);
end;
xlim([0 totalLen]);
axis tight;
ylabel(sprintf('time (%Gsec)',timebase));
xlabel('X (mm)');
xlim([0 totalLen]);
title(['exp-' num2str(dirname(1:8)) '   event-' num2str(eventNum)]);
%title(sprintf('%s, ev %d',dirname,eventNum));
colorbar
