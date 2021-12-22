function exp_details=mixVK_expDetailsRead(exp_dir)


%------------Read exp_details.txt
path_exp_details=[exp_dir '\exp_details.txt'];
fid=fopen(path_exp_details);
exp_details.factor_Master=str2num(fgetl(fid));
exp_details.factor_Slave=str2num(fgetl(fid));
exp_details.offset_Master=str2num(fgetl(fid));
exp_details.offset_Slave=str2num(fgetl(fid));
exp_details.sg_Master_order=str2num(fgetl(fid));
exp_details.sg_Slave_order=str2num(fgetl(fid));
exp_details.rosette_number=str2num(fgetl(fid));

exp_details.sg_angle=(str2num(fgetl(fid)))*pi/180;
exp_details.sg_angle=exp_details.sg_angle(exp_details.rosette_number);
exp_details.x_sg=str2num(fgetl(fid));
exp_details.x_sg=exp_details.x_sg(exp_details.rosette_number);
exp_details.sg_height=str2num(fgetl(fid));
exp_details.sg_height=exp_details.sg_height(exp_details.rosette_number);

tmp=str2num(fgetl(fid));
exp_details.UpperBlockLength=tmp(1);
exp_details.UpperBlockWidth=tmp(2);
exp_details.triggerDelay=str2num(fgetl(fid));
exp_details.MasterVref=str2num(fgetl(fid));
exp_details.MasterCurrent=str2num(fgetl(fid));
exp_details.MasterResistances=str2num(fgetl(fid));
angle(1,:)=(str2num(fgetl(fid)))*pi/180;
angle(1,:)=angle(1,exp_details.rosette_number);
angle(3,:)=(str2num(fgetl(fid)))*pi/180;
angle(3,:)=angle(3,exp_details.rosette_number);
angle(2,:)=exp_details.sg_angle;
exp_details.angle=angle;

ch_names=textscan(fid,'%s',32);
exp_details.ch_Master_names=ch_names{1};
ch_names=textscan(fid,'%s',32);
exp_details.ch_Slave_names=ch_names{1};
fclose(fid);

%-----sort by x_sg location
[exp_details.x_sg,exp_details.axis_sort_index]=sort(exp_details.x_sg);
exp_details.sg_angle=exp_details.sg_angle(exp_details.axis_sort_index);
exp_details.sg_height=exp_details.sg_height(exp_details.axis_sort_index);




