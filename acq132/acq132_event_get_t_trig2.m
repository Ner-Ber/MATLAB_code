function [M_e M_t S_e S_t]=acq132_event_get_t_trig2(exp_dir)
%M_e, S_e - the evevnt number for the master and slave respectivly. 
%M_t, S_t - t trigger for the master and slave cards, respectivly.
current_dir=pwd;
path_Master=[current_dir '\' exp_dir '\acq132_093\multivent'];
path_Slave=[current_dir '\' exp_dir '\acq132_094\multivent'];

event_dir_Master=dir([path_Master '\*.COOKED']);
event_dir_Master={event_dir_Master.name};
event_dir_Slave=dir([path_Slave '\*.COOKED']);
event_dir_Slave={event_dir_Slave.name};

%--------Master
event_dir=event_dir_Master;
path_partial=path_Master;

for k=1:length(event_dir)
    
    M_e(k)=str2num(event_dir{k}(1:2));
    path=[path_partial '\' event_dir{k}];
    M_t(k)=get_t_trig(path);
    
end

%--------Slave
event_dir=event_dir_Slave;
path_partial=path_Slave;

for k=1:length(event_dir)
    
    S_e(k)=str2num(event_dir{k}(1:2));
    path=[path_partial '\' event_dir{k}];
    S_t(k)=get_t_trig(path);
    
end

%-------Compare Master and Slave




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

