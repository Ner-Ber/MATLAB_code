function [acqSOut]=synchronizeSlowPhAcq(acqS,tVec)
%The function gets acq struct and slower t vector and returns the acq
%struct at the time of the tVec after appropriate smoothing

%the function reads acq and Phantom and returns A,acq at the same time and position
% exp=;
% Tstart=;
% Tend=;
% acqS=acq132_stream_load_mat([exp '\stream'],Tstart,'min',Tend,'Syy','Uyy','Sxy','x_sg');
% phS=phantomGetLines(exp,0,Tstart,0.5,Tend,1,9,'all',8);
% A=get_A_at_x(phS,acqS.x_sg,11);


for j=1:length(tVec)
    [~,index(j)]=min(abs(acqS.t-tVec(j)));
end

smt=floor(mean(diff(tVec(1:20)))/mean(diff(acqS.t(1:20))));%usually acq sample rate faster then camera
acqSOut=acqS;
acqSOut.t=acqS.t(index,:);

field_names=fieldnames(acqS);
for j=1:length(field_names)
    if ~(strcmp(field_names{j},'t')) && ~(strcmp(field_names{j},'exp')) && ~(strcmp(field_names{j},'x_sg')) && ~(strcmp(field_names{j},'Date') )&& ~(strcmp(field_names{j},'Path'))
    acqS.(field_names{j})=my_smooth(acqS.(field_names{j}),smt);
    acqSOut.(field_names{j})=acqS.(field_names{j})(index,:);
    end
end

