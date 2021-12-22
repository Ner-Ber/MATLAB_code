function [Lines Included] = phantomReadLines_BigPic(exp_dir,eventNum,startl,interval,endl,smt,includeFlag,lineNum)
%eventNum=0 is slow acquisition.
%lineNum -> the vector of numbers of rows from the image.
%includeFlag -> 'all' for include all.
%returns [Lines Included]

path= [exp_dir '\PhBig'];
%--------Read Phantom Meta
PhMeta = phantomReadMeta(path);

if (nargin<8 || strcmp(lineNum,'all')) %the full image is taken
    lineNum=(1:PhMeta.ImageHeight);
end

if(nargin<7 )
    includeFlag='';
end

%------
f = fopen([path '\' num2str(eventNum) '.bin'],'r');
imSize=PhMeta.ImageHeight*PhMeta.ImageWidth;

Lines=zeros(floor((endl-startl)/interval)+1,PhMeta.ImageWidth);

%------ create include mat. for '...\cal' take from last line. Because at
%the begining no Fn -> no light .Careful if there is any strong tilt

if strcmp('cal',exp_dir(end-2:end))
    Includel=endl;
else
    Includel=1;
end

fseek(f,(Includel-1)*imSize*2,'bof'); %2 for 16 bit images
Im=fread(f,imSize,'uint16');
Imr=reshape(Im',PhMeta.ImageWidth,PhMeta.ImageHeight)';
Imr=Imr(lineNum,:);

if(~strcmp(includeFlag,'all'))
Included=Imr;
Included(:,:)=0;
Included((Imr <2^12-1)&(Imr>100))=1;
else
Included=Imr;
Included(:,:)=1;
end


%------- Read lines
fseek(f,(startl-1)*imSize*2,'bof'); %2 for 16 bit images
for ind=startl:interval:endl
    
    Im=fread(f,imSize,'uint16');
    Imr=reshape(Im',PhMeta.ImageWidth,PhMeta.ImageHeight)';
    Imr=Imr(lineNum,:);
    if (smt>1)
        %Ims.Lines(floor((ind-(startl))/interval)+1,:)=smooth(mean(Imr,1),smt);
        Lines(floor((ind-(startl))/interval)+1,:)=smooth(phantomMeanLinesNoOverload(Imr,Included),smt);
    else
        %Ims.Lines(floor((ind-(startl))/interval)+1,:)=mean(Imr,1);
        Lines(floor((ind-(startl))/interval)+1,:)=phantomMeanLinesNoOverload(Imr,Included);
    end
    
    %     %we already read one line, only jump interval-1...
    if interval>1
        fseek(f,(interval-1)*imSize*2,'cof');%2 for 16 bit images
    end
    
end

fclose(f);


