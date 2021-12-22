function scan=scanplot(prefix,lines)
%read scan from disk
[scan bufferLen]=scanread(prefix);
%calculate shift via least mean squares on first two lines
shift=calibrate_scan_lines(scan,bufferLen);

%don't overdo it
lines=min(lines,length(scan));
%don't plot calibration lines
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
    plot(pos1,depth1);
  hold all
end
hold off