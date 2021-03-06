function vdsPlotSlowSliceXLim(dirname,startl,interval,endl,smt)
%oded used those parameters as function arguments
totalLen=150;
normaLine=startl;

dirname=[dirname '\vds'];

%check max left and right limits
[leftLim rightLim]=vdsCalcLims(dirname,1);
cd ..
eval(sprintf('cd %s',dirname))
ims=vdsReadSlowImsSlice(startl,interval,endl,smt);
vdst=vdsReadTimeStamps;
% if ims(1).Im == -1
%     cd ..
%     return
% end
%already smoothed each line
lines=[ims.Line]';
imsLen=length(ims);
if smt>1
    %     x=1:1:length(lines(1,:));
    %     x=totalLen-((1:1:length(lines(1,:)))-leftLim)./(rightLim-leftLim)
    %     *totalLen; %Ilya's system is mirror
    x=((1:1:length(lines(1,:)))-leftLim)./(rightLim-leftLim)*totalLen; % my system
    X=repmat(x,length(ims),1);
    t=vdst.Slow(startl:interval:startl+(imsLen-1)*interval);
    T=repmat(t',length(x),1)';
    %not really first line anymore
    firstLine=repmat(lines(normaLine,:),length(ims),1);
    %mesh(lines./firstLine);view(0,90);
    mesh(X,T,lines./firstLine,lines./firstLine);%view(0,90);
    % mesh(lines./firstLine);view(0,90);
else
    %     x=1:1:length(lines(1,:));
    %     x=totalLen-((1:1:length(lines(1,:)))-leftLim)./(rightLim-leftLim)
    %     *totalLen; %Ilya's system is mirror
    x=((1:1:length(lines(1,:)))-leftLim)./(rightLim-leftLim)*totalLen; % my system
    X=repmat(x,length(ims),1);
    t=vdst.Slow(startl:interval:startl+(imsLen-1)*interval);
    T=repmat(t',length(x),1)';
    %not really first line anymore
    %     firstLine=repmat(smooth(lines(normaLine,:),10)',length(ims),1);
    %     mesh(lines);view(0,90);
    mesh(X,T,lines);view(0,90);
end;
%ylabel('time (sec)');
xlim([0 totalLen]);
ylabel('time (sec)');
xlabel('x (mm)');
colorbar
view(0,90);

cd ..
cd ..