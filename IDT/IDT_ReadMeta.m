function outdata=IDT_ReadMeta(path)
%full path of the directories
if nargin<1
    path='';
end
IDT_dir = dir([path,'\*xsv']);

%--- read cih file from last event:
outdata = Read_xsv_file([path '\' IDT_dir(end).name]);

%--- check if frameRate makes sense (may not make of external clock, then
%the default will be fps=200,000)
if outdata.FrameRate==1
    outdata.FrameRate=200000;
end

%--- get number of events
FieldValues = struct2cell(outdata);
outdata.NumEvents = nnz(cell2mat(FieldValues(~cellfun(@isempty,regexp(fieldnames(outdata),'BROC_Time_')))));

outdata.PostIms = outdata.BROCLenght - outdata.PerTrigFrms;
outdata.FrameT=(1/outdata.FrameRate)*1e6; %[musec]

Start_frame=outdata.StartFrame;
Frame_rate=outdata.FrameRate;
trig_delay=7e-3;                       % trig delay is 7 msec
event_frame=round(Start_frame-trig_delay*Frame_rate);
outdata.EventFrame=event_frame;


