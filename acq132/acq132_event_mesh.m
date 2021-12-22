function acq132_event_mesh(exp_dir,event,objName,Tstart,Tend,smt,normalize,baseFig)
%event 0 is slow 
%pre,post are in msec in refernce to the trig.
%normalize=1 -> FirstLine is subtructed:               D=D-FirstLine;
%normalize=2 -> The data normolized by FirstLine:      D=D./FirstLine;


if (nargin<8)
    baseFig=gcf;
end

if(event==0)
    acq=acq132_stream_load_mat(exp_dir,Tstart,0.5,Tend,'x_sg',objName);
    obj=eval(sprintf(['acq.' objName]));
else
    acq=acq132_event_get_data(exp_dir,event,Tstart,Tend,smt,'x_sg',objName);
    obj=eval(sprintf(['acq.' objName]));
end
set(0,'DefaultFigureWindowStyle','docked')
%------------------meshes

figure(baseFig);
my_mesh(acq.x_sg([1:5 7:15]),acq.t,obj(:,[1:5 7:15]),normalize); %sg number 6 (150mm block)is not wroking!!!!!!!!!!!!!!! 
title([exp_dir(1:8) ' event-' num2str(event) ' Ttrigger-' num2str(acq.t_trig) ' ' objName]);
xlim([0 150]);
