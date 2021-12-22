function outdata=readVdsMeta
txtData=dlmread(sprintf('vds_meta.txt'));
outdata.NumEvents=txtData(1);
outdata.EventNumBuffers=txtData(2);
outdata.BufferNumImages=txtData(3);
outdata.ImageHeight=txtData(4);
outdata.ImageWidth=txtData(5);
% CMC 1300 takes 2\musec for each line plus 1 blanking line
outdata.FrameT=(outdata.ImageHeight+1)*2E-6;
outdata.FrameRate=1/outdata.FrameT;
outdata.PostIms=txtData(7);
