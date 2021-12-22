function outt=vdsReadTimeStamps
Rapid=readVdsRapidTimeStamps;
Slow=readVdsSlowTimeStamps;
%read meta data to eliminate offset 
vdsMeta = readVdsMeta;
roffset=vdsMeta.FrameT*vdsMeta.PostIms;


%aline all timestamps to the first "slow clock" tick
if (length(Slow)>10)
    p=polyfit(1:10,Slow(1:10)',1);
else
    p=polyfit(1:length(Slow),Slow,1);
end
outt.Slow=Slow-p(2);
outt.Rapid=Rapid-p(2)-roffset;
if (outt.Rapid(1)<0)
    fileTimes=dlmread('VdsRapidTimeStamps.txt');
    %first row is the frequency of the counter in counts per second
    freq=fileTimes(1);
    outt.Rapid=outt.Rapid+2^32/freq;
end
