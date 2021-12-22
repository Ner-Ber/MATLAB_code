function Ims=phantomGetLinesStd(exp_dir,eventNum,Tstart,Tinterval,Tend,timeBase,smt,includeFlag,lineNum)
%Ims=phantomGetLines(exp_dir,eventNum,Tstart,Tinterval,Tend,timeBase,smt,includeFlag,lineNum)
%eventNum=0 is slow acquisition. timeBase=1000 -> msec
%returns  Ims.Exp Ims.t_trigger  Ims.t Ims.Lines Ims.Included Ims.x
%if a.tiff exists Ims is normalized


if (nargin==5)
    timeBase=1;
    smt=0;
elseif (nargin==6)
    smt=0;
end

if (nargin<8)
    includeFlag='';
end
if (nargin<9)
    lineNum='all';
end

Ims.Exp=exp_dir;
Ims.EventNum=eventNum;

path= [exp_dir '\Ph'];
%-------trigger Time
if (eventNum~=0)
    t_trigger=dlmread([path '\eventPhTriggerTime.txt']);
    Ims.t_trigger=t_trigger(eventNum);
end

[startl,interval,endl,Ims.t]=phantomReadTime(exp_dir,eventNum,Tstart,Tinterval,Tend,timeBase); %this function is responsible for all the time logic
%[Ims.lines Ims.Included]=phantomReadLinesContrast(exp_dir,eventNum,startl,interval,endl,smt,includeFlag,lineNum);
[Ims.lines Ims.Included]=phantomReadLinesStd(exp_dir,eventNum,startl,interval,endl,smt,includeFlag,lineNum);

Ims.x=phantomBuildXaxis(exp_dir(1:8)); %if there is tilt can't calculate limits beacause not all the block is presses. -> calculate limits after pushing started
