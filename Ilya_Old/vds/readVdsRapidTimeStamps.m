function t=readVdsRapidTimeStamps
fileTimes=dlmread('VdsRapidTimeStamps.txt');
%first row is the frequency of the counter in counts per second
freq=fileTimes(1);
%get event times in seconds
t=fileTimes(2:length(fileTimes))/freq;
% dt=diff(t);
% ind=1:length(t);
% indt=ind(dt<0)+1;
% t(ind(indt:length(t)))=t(ind(indt:length(t)))+2^32/freq;
% 
% 
% if(t(1)<0)
%     t=t+2^32/freq;
% end

dt=diff(t);
ind=1:length(t);
indt=ind(dt<0)+1;
if (length(indt)<2)
  t(ind(indt:length(t)))=t(ind(indt:length(t)))+2^32/freq;
else
    for j=1:length(indt)
        t(ind(indt(j)):length(t))=t(ind(indt(j)):length(t))+2^32/freq;
    end
end 