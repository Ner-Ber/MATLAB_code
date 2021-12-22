function [leftLim rightLim]=phantomCalcLims(exp_dir,Tstart,Tinterval,Tend)
  
[startl,interval,endl,Ims.t]=phantomReadTime(exp_dir,0,Tstart,Tinterval,Tend,1); %this function is responsible for all the time logic

%The function gets the slow lines 

lines1=phantomReadLines(exp_dir,0,startl,interval,endl,1,'all',2:7);


% subtruct first line for offset correction . (needed if there any background light)
% if strcmp(exp_dir(end-2:end),'cal') || strcmp(exp_dir(end-3:end),'cal\') %if path is cal directory so this is also cal_dir
%     cal_dir=exp_dir;
% elseif exist([exp_dir '\cal'],'dir')==7 % if cal directory exists offset correction taken from there, otherwise from the first event
%     cal_dir=[exp_dir '\cal'];
% else
%     cal_dir=exp_dir;
% end
% line0=phantomReadLines(cal_dir,0,1,1,1,1,'all');
% lines1=lines1-repmat(line0,length(lines1(:,1)),1);

left=zeros(length(lines1(:,1)),1)+NaN;
indLeft=1;
right=zeros(length(lines1(:,1)),1)+NaN;
indRight=1;

for i=1:length(lines1(:,1))
    
    buffer=find(lines1(i,:)>300,1,'first');
    if(any(buffer))
    left(indLeft)=buffer;
    indLeft=indLeft+1;
    end
    buffer=find(lines1(i,:)>300,1,'last');
    
    if(any(buffer))
    right(indRight)=buffer;
    indRight=indRight+1;
    end
end

rightLim=max(smooth(right,100));
leftLim=min(smooth(left,100));


