function outIms = vdsReadImsSlice(eventNum,startl,interval,endl,smt)
vdsMeta = readVdsMeta;
if eventNum>vdsMeta.NumEvents
    sprintf('Event Num too high')
    outIms(1).Ind =-1;
    cd ..
    return
end
filename=sprintf('%d.bin',eventNum-1);
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

  if ind==startl-1
            Included=Imr;
             Included(:,:)=0;
             Included(Imr <255 & Imr>50) =1;
%             Included(:,:)=1;
 end
% Included=Imr;
% Included(:,:)=1;
%    mesh(Imr');view(0,90);
%    pause

    if (smt>1)
%         Ims((ind-(startl-1))/interval+1).Line=smooth(mean(Imr,2),smt);
            Ims((ind-(startl-1))/interval+1).Line=smooth(meanLinesNoOverload(Imr,Included),smt);
    else
%         Ims((ind-(startl-1))/interval+1).Line=mean(Imr,2);
            Ims((ind-(startl-1))/interval+1).Line=meanLinesNoOverload(Imr,Included);
    end
    Ims((ind-(startl-1))/interval+1).Ind=ind;
    %we already read one line, only jump interval-1...
    if interval>1
        fseek(f,(interval-1)*imSize,'cof');
    end
end
fclose(f);
outIms=Ims;