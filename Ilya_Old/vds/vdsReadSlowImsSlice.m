function outIms = vdsReadSlowImsSlice(startl,interval,endl,smt)
vdsMeta = readVdsMeta;
filename='SlowVds.bin';

f = fopen(filename);
imSize=vdsMeta.ImageHeight*vdsMeta.ImageWidth;
fseek(f,(startl-1)*imSize,'bof');
% ind=1;
for ind=startl-1:interval:endl
  Im=fread(f,imSize,'uchar');
  if (feof(f)) 
       %lim was too far, stop reading
%        sprintf('Lim too high %d - cut at end - line %d',endl,ind)
       break;
   end
  Imr=reshape(Im,vdsMeta.ImageWidth,vdsMeta.ImageHeight);
%    mesh(Imr');view(0,90);
%    pause
    if (smt>1)
        Ims((ind-(startl-1))/interval+1).Line=smooth(mean(Imr,2),smt);
    else
        Ims((ind-(startl-1))/interval+1).Line=mean(Imr,2);
    end
    
    %we already read one line, only jump interval-1...
    if interval>1
        fseek(f,(interval-1)*imSize,'cof');
    end
end
fclose(f);
outIms=Ims;