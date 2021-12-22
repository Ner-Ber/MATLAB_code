function phantomPlotSliceXlim(exp_dir,eventNum,Tstart,Tinterval,Tend,timeBase,smt,lineNum,baseFig)
%eventNum=0 is slow acquisition. [msec] timeBase=1000.


if nargin<9
    baseFig=gcf;
end
if nargin<8
    lineNum='all';
end
set(0,'DefaultFigureWindowStyle','docked')
exp_details=expDetailsRead(exp_dir);

ims=phantomGetLines(exp_dir,eventNum,Tstart,Tinterval,Tend,timeBase,smt,'all',lineNum);

[~,index1]=min(abs(ims.x-6));
[~,index2]=min(abs(ims.x-exp_details.UpperBlockLength));

figure(baseFig);
my_mesh(ims.x(index1:index2),ims.t,ims.lines(:,(index1:index2)),2);
%imagesc(ims.x(index1:index2),ims.t,ims.lines(:,(index1:index2)));
set(gca,'ydir','normal')
% 
 caxis([0.6 1.1]);
% xlabel('X (mm)');
 xlim([0 exp_details.UpperBlockLength]);
 title([ ims.Date ' ' ims.Exp '  event - ' num2str(eventNum) '  Trigger=' num2str(ims.t_trigger)]);


%title(['exp-' num2str(dirname(1:8)) '   event-' num2str(eventNum)]);
%title(sprintf('%s, ev %d',dirname,eventNum));
%timebase=1; %[msec]
%ylabel(sprintf('time (%Gsec)',timebase));
%colorbar
