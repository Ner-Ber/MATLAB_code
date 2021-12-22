function vdsPlotEvent_x_sg(exp_dir,eventNum,startl,endl,smt,BaseFig)
%The function plotes the contact area at sg locations

if nargin<6
    BaseFig=1;
end

normaLine=1;
interval=1;
totalLen=150;

[t_vds,x,line]= vds_event_get_data(exp_dir,eventNum,startl,endl,smt);

%------------Read exp_details.txt
path_exp_details=[exp_dir '\exp_details.txt'];
my_format=dlmread(path_exp_details, ' ' ,[0,0,5,31]);
x_sg=150-nonzeros(my_format(6,:))'; %notice the positive x axis on the exp_details is in the pushing direction. we work in the camera frame (opposite direction)
x_sg=sort(x_sg);

pix=zeros(1,length(x_sg));

for sg_num=1:15
    [~,pix(sg_num)]=min(abs(x-x_sg(sg_num)));
end

FirstLine=repmat(mean(line(1:5,pix)),length(t_vds),1);
figure(BaseFig);
hold on
plot(t_vds,line(:,pix)./FirstLine-1,'o-')

