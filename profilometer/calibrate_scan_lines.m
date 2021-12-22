function shift = calibrate_scan_lines(scan,bufferLen)
%return shift
%indata=indata - mean(mean(indata,1));
% shift = fminbnd(@(p) crosscorr_lines(p,indata),-max,-1);

max=bufferLen;

for i=-max:max
    msq(i+max+1)=crosscorr_lines(i,scan);
end
msq(msq==0)=nan;
[m shift]=min(msq);
shift=shift-1-max;
