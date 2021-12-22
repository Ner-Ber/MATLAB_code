clear
dirName='12-50-17';
event=5;


% ------slowStream
% phs=PhStreamRead([dirName '\Ph'],1,1,'max',0);
% phs.Line= sum(phs.Lines,2);
% phs.Line=phs.Line/max(phs.Line);
% 
% [acqs acqsT]=acq132_stream_read([pwd '\' dirName '\acq132_094\stream']);
% acqs=acqs(:,[16 31]);
% acqs=acqs-repmat(mean(acqs(100:150,:)),length(acqs(:,1)),1);
% acqs(:,1)=abs(acqs(:,1));
% acqs=acqs./repmat(max(acqs),length(acqs(:,1)),1);
% 
% figure(1); 
% plot(acqsT,acqs,'.-');
%  hold on
%  plot(phs.t+phs.t(1)/2,phs.Line,'ro-');
%  hold off
 %-------

%phEventPlotSliceXlim(dir,event,'max',50,'max',0);

eval(sprintf('cd %s',[dirName '\Ph']));
phe=PhEventReadImsSlice(event,'max',1,'max',1);
cd ..
cd ..
phe.Line=sum(phe.Lines,2);
phe.Line=phe.Line/max(phe.Line);


[acqe,acqeT,acqeTtrig]=acq132_event_read([pwd '\' dirName '\acq132_094\multivent'],event);
acqe=acqe(:,[16 31]);
acqe(:,1)=acqe(:,1)/max(acqe(:,1));
acqe(:,2)=acqe(:,2);max(acqe(:,2));


 figure(2);
 plot(acqeT,acqe,'.-');
 hold all;
 for l=400:50:450
 plot(phe.t,phe.Lines(:,l)/max(phe.Lines(:,l)),'o-');
 
 end
 hold off
 
