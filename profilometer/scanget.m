function ss=scanget(prefix,lines,xStart,xEnd)

%if xStart and xEnd are specified one x axis is created and all data is
%interpolated to that axis;

%note: lines includes the 2 calibration lines. for 6 lines enter lines=8

%read scan from disk
[scan bufferLen]=scanread(prefix);
%calculate shift via least mean squares on first two lines
shift=calibrate_scan_lines(scan,bufferLen);

%don't overdo it
lines=min(lines,length(scan));

for i=1:lines
    if(shift>0)
        pos1=scan(i).pos(shift:end);
        depth1=scan(i).depth(1:end-shift+1);
    elseif(shift<0)
        pos1=scan(i).pos(1:end+shift+1);
        depth1=scan(i).depth(-shift:end);
    else
        pos1=scan(i).pos;
        depth1=scan(i).depth;
    end
    
    
    %sometimes the stage doesnt move -> neglect thos points
    diffPos=diff(pos1);
    if (mean(diffPos)>0)
        pos1=pos1(diffPos>0);
        depth1=depth1(diffPos>0);
    else
        pos1=pos1(diffPos<0);
        depth1=depth1(diffPos<0);
    end
    
    s{i}.pos=pos1(1:end-1);
    s{i}.depth=depth1(1:end-1);
    sLength(i)=length(s{i}.pos);
end

%don't plot calibration lines
minsLength=min(sLength(3:end));
for i=3:lines
        ss.pos(:,i-2)=s{i}.pos(1:minsLength)';
        ss.depth(:,i-2)=s{i}.depth(1:minsLength)';
end

%if boundries are specified one x axis is created and all data is
%interpolated to that axis;
if (nargin>2) 
    
    x=logical((ss.pos>xStart-10).*(ss.pos<xEnd+10));
    ssCut.pos=(xStart:xEnd)';
    ssCut.depth=zeros(length(ssCut.pos),size(ss.depth,2));
    
    for i=1:length(ss.pos(1,:))    
        ssCut.depth(:,i)=interp1(ss.pos(x(:,i),i),ss.depth(x(:,i),i),ssCut.pos);
    end
    
    ss=ssCut;
end
