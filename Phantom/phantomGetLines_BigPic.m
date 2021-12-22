function Ims=phantomGetLines_BigPic(exp_dir,eventNum,Tstart,Tinterval,Tend,timeBase,smt,includeFlag,lineNum)
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

path=pwd;
lastSlashPosition = find(path == '\', 1, 'last');
Ims.Date = path(lastSlashPosition+1:end);

Ims.Exp=exp_dir;
Ims.EventNum=eventNum;

path= [exp_dir '\PhBig'];
%-------trigger Time
if (eventNum~=0)
    t_trigger=dlmread([path '\eventPhTriggerTime.txt']);
    Ims.t_trigger=t_trigger(eventNum);
end

[startl,interval,endl,Ims.t]=phantomReadTime_BigPic(exp_dir,eventNum,Tstart,Tinterval,Tend,timeBase); %this function is responsible for all the time logic
[Ims.lines Ims.Included]=phantomReadLines_BigPic(exp_dir,eventNum,startl,interval,endl,smt,includeFlag,lineNum);

exp_details=expDetailsRead(exp_dir);
L=0.004; exp_details.UpperBlockLength;
% path=[exp_dir '\Ph\a.tif'];
% if(exist(path))
%     a=double(imread(path));%/2^16*2^12;
%     a=mean(a,1);
%     %line0=smooth(a,7)';
%     x=1:length(a); 
%     p=polyfit(x,a,3);
%     line0=polyval(p,x);
%     Ims.lines=Ims.lines./repmat(line0,length(Ims.lines(:,1)),1);
% end

%subtruct first line for offset correction . (needed if there any background light)
% 
% if strcmp(exp_dir(end-2:end),'cal') || strcmp(exp_dir(end-3:end),'cal\') %if path is cal directory so this is also cal_dir
%     cal_dir=exp_dir;
% elseif exist([exp_dir '\cal'],'dir')==7 % if cal directory exists offset correction taken from there, otherwise from the first event
%     cal_dir=[exp_dir '\cal'];
% else
%     cal_dir=exp_dir;
% end
% [line0 ~]=phantomReadLines(cal_dir,1,1,1,1,1,includeFlag,lineNum);
% 
% Ims.lines=Ims.lines-repmat(line0,length(Ims.lines(:,1)),1);


%if the contact are starts at pix1 and end at pix2
   pix1=1;%9;%27;
   pix2=size(Ims.lines,2);%1263;%1267;
    %pix1=5;
    %pix2=1255;
%   Ims.x=(-pix1+1:1280-pix1)/(pix2-pix1-1)*L; %When there is no slow
% % 
%     pix1=190;
%     pix2=1180;
  Ims.x=(-pix1+1:size(Ims.lines,2)-pix1)/(pix2-pix1-1)*L; %When there is no slow
%Ims.x=phantomBuildXaxis(exp_dir(1:8)); %if there is tilt can't calculate limits beacause not all the block is presses. -> calculate limits after pushing started
