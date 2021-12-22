function exp_details=expDetailsRead_till_2012_10_24(exp_dir)
%------------Read exp_details.txt

%2 Very important and confusing notes:
%1.Because of the routiong betwean the rossetes channels and amplifing and
%Peter's channel numbering I'm using sg_093_order and sg_094_order to
%arrange the channels from pushing side to the stoppper .
%2.Thus in exp_details.txt the positive x directionis from the pushing side to the stoppper direction. Ussually we
%use opposite direction. I converted the x_sg units and angle sign and the
%order. Use exp_details.axis_sort_index to sort the strains in acq132_convert_raw_to_data.

LastLineNum=10;
path_exp_details=[exp_dir '\exp_details_till_2012_10_24.txt'];
my_format=dlmread(path_exp_details, ' ' ,[0,0,LastLineNum,31]);
exp_details.factor_093=my_format(1,:);
exp_details.factor_094=my_format(2,:);
exp_details.offset_093=my_format(10,:);
exp_details.offset_094=my_format(11,:);
exp_details.sg_093_order=nonzeros(my_format(3,:))';
exp_details.sg_094_order=nonzeros(my_format(4,:))';
exp_details.UpperBlockLength=my_format(8,1);%x dimension
exp_details.UpperBlockWidth=my_format(8,2);
exp_details.triggerDelay=nonzeros(my_format(9,:));

exp_details.x_sg=exp_details.UpperBlockLength-nonzeros(my_format(6,:))'; %notice the positive x axis on the exp_details is in the pushing direction.Later we work in the camera frame (opposite direction)
[exp_details.x_sg,exp_details.axis_sort_index]=sort(exp_details.x_sg);

exp_details.sg_angle=-(nonzeros(my_format(5,:)))'*pi/180;%the "-" is needed to change the x axis direction  [radians].
exp_details.sg_angle=exp_details.sg_angle(exp_details.axis_sort_index);
exp_details.sg_height=nonzeros(my_format(7,:))';%notice the positive x axis on the exp_details is in the pushing direction.Later we work in the camera frame (opposite direction)
exp_details.sg_height=exp_details.sg_height(exp_details.axis_sort_index);

%read Channel names---allways keep at the end
fid=fopen(path_exp_details);
ch_names=textscan(fid,'%s',32,'HeaderLines',LastLineNum+2);
exp_details.ch_Master_names=ch_names{1};
ch_names=textscan(fid,'%s',32);
exp_details.ch_Slave_names=ch_names{1};
fclose(fid);
