function plot_Sxy_Ds_A(exp_dir,event,sg_num,ch_094_ex_num,baseFig)
%plot_Sxy_Ds_A(exp_dir,event,sg_num,ch_094_ex_num,baseFig)
%The function calculates the Stresses,displacement and A at perticular sg.
%And normolize the data. 

if nargin<5
    baseFig=gcf;
end


e=acq132_event_get_data(exp_dir,event,'max','max',1,'Sxx','Sxy','Syy','ch_094_ex','x_sg');
f=calc_Fx_stress_splined(e);
index= (f.x==e.x_sg(sg_num));
Fx=f.Fx(:,index);
clear f;

% Ods=e.ch_094_ex(:,ch_094_ex_num); %displacment in lab frame 
% ds=e.ch_094_ex(:,3);
% slip=Ods;
%slip=Ds-(F-mean(F(1:20)))/sensF; %slip relative to the stage

%---Cut single pix v.s t

if  exist([exp_dir '\vds'],'dir')==7
vds=vds_event_get_data(exp_dir,event,1,1000,3);
[~,pix]=min(abs(vds.x-e.x_sg(sg_num)));
line=subtruct_norm(vds.lines(:,pix),1); %normilize the line;
vds=rmfield(vds,'x');
vds=rmfield(vds,'lines');
end

if  exist([exp_dir '\Ph'],'dir')==7
vds=vds_event_get_data(exp_dir,event,1,1000,3);
[~,pix]=min(abs(vds.x-e.x_sg(sg_num)));
line=subtruct_norm(vds.lines(:,pix),1); %normilize the line;
vds=rmfield(vds,'x');
vds=rmfield(vds,'lines');
end

%-----plot

% figure(baseFig)
% plotyy(e.t,e.Sxy(:,sg_num) ,e.t,[slip ds])
% baseFig=baseFig+1;
% my_xlim([e.t(1) e.t(end)])
%title(['ch094ex' num2str(ch_094_ex_num) ',num2str(e.x(sg) ])

figure(baseFig);
plotyy(e.t,e.Sxy(:,sg_num),vds.t,line)
baseFig=baseFig+1;
my_xlim([e.t(1) e.t(end)]) 

figure(baseFig);
plotyy(e.t,e.Sxy(:,sg_num),e.t,Fx)
my_xlim([e.t(1) e.t(end)]) 
