function t_trig=acq132_event_get_t_trig(exp_dir)
% old version - in some cases doesn't work properly. Use
% acq132_event_get_t_trig2

current_dir=pwd;
path_Master=[current_dir '\' exp_dir '\acq132_093\multivent'];
path_Slave=[current_dir '\' exp_dir '\acq132_094\multivent'];

event_dir_Master=dir([path_Master '\*.COOKED']);
event_dir_Slave=dir([path_Slave '\*.COOKED']);

%--Sometime one of the cards misses an event
if length(event_dir_Master)==length(event_dir_Slave)
    event_dir={event_dir_Master.name};
    path_partial=path_Master;
elseif length(event_dir_Master)<length(event_dir_Slave)
    event_dir={event_dir_Slave.name};
    path_partial=path_Slave;
    display('Master card missed an event')
elseif length(event_dir_Master)>length(event_dir_Slave)
    event_dir={event_dir_Master.name};
    path_partial=path_Master;
    display('Slave card missed an event')
end

for k=1:length(event_dir)
    
    path=[path_partial '\' event_dir{k}];
    if (str2num(event_dir{k}(1:2))-k ~= 0)
        display(['Both cards missed event ' num2str(k)])
    end
    
    t_trig(k)=get_t_trig(path);
    
end

function t_trig=get_t_trig(path)

%-----------find pre
newData1 = importdata([path '\format'],'\t',42);
str=newData1.textdata{3,1};
s1=strfind(str,'--pre');
pre=sscanf(str(s1:end),'--pre %d');

%----------------Time base read--------
file=dir([path '\*.TBD']);
fid= fopen([path '\' file.name],'r');
t=fread(fid, 'float64');
t_trig=t(pre-2)+(t(pre-1)-t(pre-2))/2;%[sec] Strange!! but was checked. would expect t_trig=t(pre)+(t(pre+1)-t(pre))/2;
%t_trig=t(pre-1)+(t(pre)-t(pre-1))/2;%[sec]
%t_trig=t(pre)+(t(pre+1)-t(pre))/2;%[sec]
fclose(fid);