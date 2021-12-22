function Ims=phantomGetIms(exp_dir,eventNum,Tstart,Tinterval,Tend,timeBase,smt,lineNum)
%Ims=phantomGetIms(exp_dir,eventNum,Tstart,Tinterval,Tend,timeBase,smt,lineNum)
%eventNum=0 is slow acquisition. timeBase=1000 -> msec
%lineNum -> the vector of numbers of rows from the image.
%returns  Ims.Exp Ims.t_trigger  Ims.t Ims.ims Ims.Included Ims.x
%if a.tiff exists Ims is normalized

if (nargin==5)
    timeBase=1;    smt=1;
elseif (nargin==6)
    smt=1;
end
if (nargin<8)
    lineNum='all';
end

current_dir=pwd;
lastSlashPosition = find(current_dir == '\', 1, 'last');
Ims.Date=current_dir(lastSlashPosition+1:end);

Ims.Exp=exp_dir;
Ims.EventNum=eventNum;

path= [exp_dir '\Ph'];
%-------trigger Time
if (eventNum~=0)
    t_trigger=dlmread([path '\eventPhTriggerTime.txt']);
    Ims.t_trigger=t_trigger(eventNum);
end

[startl,interval,endl,Ims.t]=phantomReadTime(exp_dir,eventNum,Tstart,Tinterval,Tend,timeBase); %this function is responsible for all the time logic
Ims.ims=phantomReadIms(exp_dir,eventNum,startl,interval,endl,smt,lineNum);

% path=[exp_dir '\Ph\a.tif'];
% if(exist(path))
%     a=double(imread(path));
%     a=mean(a(1:8,:),1);
%     %line0=smooth(a,7)';
%     x=1:length(a); 
%     p=polyfit(x,a,4);
%     line0=polyval(p,x);
%     Ims.ims=Ims.ims./repmat(line0,[size(Ims.ims,1),1,size(Ims.ims,3)]);
% end

Ims.x=phantomBuildXaxis(exp_dir(1:8)); %if there is tilt can't calculate limits beacause not all the block is presses. -> calculate limits after pushing started



