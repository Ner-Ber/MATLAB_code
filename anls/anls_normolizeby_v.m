set(0,'DefaultFigureWindowStyle','docked')

clear
clc
pre=7.6;
post=-5.01;
baseFig=1;
exp_dir='17-26-41';
event=15;
sg_num=4;

smt=7;

[t_trig,t_acq,x_sg,~,~,~,Sxx,Syy,Sxy,N,F]=acq132_event_get_data(exp_dir,event,pre,post,smt,{'Sxy','Sxx','Syy'});

[sg t_front pk_front v_front]=acq132_event_calc_v(exp_dir,event,pre,post,smt);
FirstLine=repmat(mean(Sxy(1:20,:)),length(t_acq),1);
D=Sxy./FirstLine-1;
figure;

for j=1:length(sg)
d=(t_acq-t_front(j))*v_front(j);
D_norm=D(:,sg(j))/pk_front(sg(j));
plot(d,D_norm,'.-');
hold all
end