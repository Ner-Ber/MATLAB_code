function outdata=photronReadMeta(path)
%full path of the directories
if nargin<1
    path='';
end
Photron_dir = dir(path);

%--- read cih file from last event:
[outdata] = Read_cih_file([path '\' Photron_dir(end).name '\' Photron_dir(end).name]);
outdata.NumEvents=size(Photron_dir,1)-2;
% outdata.PostIms = outdata.NumIms - outdata.StartFrame;
outdata.PostIms = outdata.TotalFrame + outdata.StartFrame;
outdata.FrameT=outdata.ShutterSpeed_s*1e6; %[musec]
outdata.FrameRate = outdata.RecordRate_fps;


%% determine event frame, starting frame and ending frame
% Start_frame=0;
% Frame_rate=PhMeta.FrameRate;
% trig_delay=7e-3;                       % trig delay is 7 msec
% event_frame=round(Start_frame-trig_delay*Frame_rate)+abs(PhMeta.StartFrame);
% PhMeta.EventFrame=event_frame;

Start_frame=0;
Frame_rate=outdata.RecordRate_fps;
trig_delay=7e-3;                       % trig delay is 7 msec
event_frame=round(Start_frame-trig_delay*Frame_rate)+abs(outdata.StartFrame);
outdata.EventFrame=event_frame;

