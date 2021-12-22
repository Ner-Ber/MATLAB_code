function outIms = vdsReadImsSlow(interval)
vdsMeta = readVdsMeta;

filename='SlowVds.bin';
f = fopen(filename);
imSize=vdsMeta.ImageHeight*vdsMeta.ImageWidth;
i=1;

while feof(f) ==0
    Im=fread(f,imSize,'uchar');
    if ~isempty(Im)
%         Ims(i).Im=reshape(Im,vdsMeta.ImageWidth,vdsMeta.ImageHeight);
        Imr=reshape(Im,vdsMeta.ImageWidth,vdsMeta.ImageHeight);
        Ims(i,:).Line1=mean(Imr(:,1:vdsMeta.ImageHeight),2);
        if interval>1
            fseek(f,(interval-1)*imSize,'cof');%jump to next after interval
        end
    end
    i=i+1;
end

outIms=Ims;