function [scan bufferLend] = scanread(prefix)
% read Buffers file
f = fopen(sprintf('%s_Buffers.bin',prefix));
i=1;
while (~feof(f))
    % read the row number and switch from zero based C to one base Matlab
    buffer=fread(f,1,'int32');
    if(~isempty(buffer))
        buffers(i)=buffer;
    end
    i=i+1;
end
fclose(f);


% read Depth file
fd = fopen(sprintf('%s_Depth.bin',prefix));
bufferLend=fread(fd,1,'int32');

% read Position file
fp = fopen(sprintf('%s_Position.bin',prefix));
bufferLenp=fread(fp,1,'int32');

for i=1:length(buffers)
    % read the row
    depth=fread(fd,bufferLend,'float32');
    pos=fread(fp,bufferLenp,'float32');

        sim(i).depth=depth;
        sim(i).pos=pos;
        rowBuffers(i)=buffers(i);
end
fclose(fd);
fclose(fp);

%now loop over the rows and split the read buffers into the scan output
%struct

%start with calibration back and forth rows -2 and -1
for i=-2:-1
        extra=-mod(i+1,2)*0.5; 
        depth=[sim(rowBuffers==i).depth];
        scan(i+3).depth=reshape(depth,size(depth,1)*size(depth,2),1);
        pos=[sim(rowBuffers==i).pos];
        scan(i+3).pos=reshape(pos,size(pos,1)*size(pos,2),1)+extra;
end

%than add the real buffers
for i=1:max(buffers)
        extra=-mod(i,2)*0.5;
        depth=[sim(rowBuffers==i).depth];
        scan(i+2).depth=reshape(depth,size(depth,1)*size(depth,2),1);
        pos=[sim(rowBuffers==i).pos];
        scan(i+2).pos=reshape(pos,size(pos,1)*size(pos,2),1)+extra;
end

