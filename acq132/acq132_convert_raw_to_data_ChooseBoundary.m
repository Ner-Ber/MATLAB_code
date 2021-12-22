function acq132=acq132_convert_raw_to_data_ChooseBoundary(exp_dir,ch_Master,ch_Slave,varargin)
% acq132=acq132_convert_raw_to_data_ChooseBoundary(exp_dir,ch_Master,ch_Slave,PlaneFlag)
% PlaneFlag is optional. default='PlaneStrain'
%
% uses "exp_details.txt" to convert units , correct offset and s.g rotation, calculate stress and strain.
%It makes the final data
%exp_dir is needed to use event1 for offset substruction

PlaneFlag = setDefaults4function(varargin,'PlaneStrain');

%------------Read exp_details.txt
exp_details=expDetailsRead(exp_dir);
factor_Master=exp_details.factor_Master;
offset_Master=exp_details.offset_Master;
offset_Slave=exp_details.offset_Slave;
factor_Slave=exp_details.factor_Slave;
sg_Master_order=exp_details.sg_Master_order;
sg_Slave_order=exp_details.sg_Slave_order;
sg_angle=exp_details.sg_angle;
x_sg=exp_details.x_sg;
y_sg=exp_details.sg_height;
ch_Master_names=exp_details.ch_Master_names;
ch_Slave_names=exp_details.ch_Slave_names;
%---------------offset correction
dirCell=path2cell(exp_dir);
calFlag=strcmp(dirCell,'cal');
cal_dir='';
if (sum(calFlag)) %if path is cal directory so this is also cal_dir
    for j=1:find(calFlag==1) cal_dir=[cal_dir '\' dirCell{j}]; end
else
    cal_dir=[dirCell{1} '\cal'];
end

% if cal directory exists offset correction taken from there, otherwise from the first event
if  exist(cal_dir,'dir')==7
    path_Master=[cal_dir '\acq132_093\multivent'];
    path_Slave=[cal_dir '\acq132_094\multivent'];
else
    path_Master=[exp_dir '\acq132_093\multivent'];
    path_Slave=[exp_dir '\acq132_094\multivent'];
    cal_dir=(exp_dir);
end

%--offset for Master
ch_offset=acq132_event_read(path_Master,1);%using event1 for offset correction
ch_offset=mean(ch_offset);
ch_number_for_offset = find(offset_Master==0);
ch_Master(:,ch_number_for_offset)=ch_Master(:,ch_number_for_offset)-repmat(ch_offset(ch_number_for_offset),length(ch_Master(:,1)),1);%offset correction// comment this line for disable
ch_number_for_offset = find(offset_Master~=0);
ch_Master(:,ch_number_for_offset)=ch_Master(:,ch_number_for_offset)-repmat(offset_Master(ch_number_for_offset),length(ch_Master(:,1)),1);%offset correction with specified values// comment this line for disable
ch_Master=ch_Master./repmat(factor_Master,length(ch_Master(:,1)),1);%unit conversion notice the using "./" instead ".*", empty channels will be "inf"

%--offset for Slave
ch_offset=acq132_event_read(path_Slave,1);%using event1 for offset correction
ch_offset=mean(ch_offset);
ch_number_for_offset = find(offset_Slave==0);
ch_Slave(:,ch_number_for_offset)=ch_Slave(:,ch_number_for_offset)-repmat(ch_offset(ch_number_for_offset),length(ch_Slave(:,1)),1);%offset correction// comment this line for disable
ch_number_for_offset = find(offset_Slave~=0);
ch_Slave(:,ch_number_for_offset)=ch_Slave(:,ch_number_for_offset)-repmat(offset_Slave(ch_number_for_offset),length(ch_Slave(:,1)),1);%offset correction with specified values// comment this line for disable
ch_Slave=ch_Slave./repmat(factor_Slave,length(ch_Slave(:,1)),1);%unit conversion notice the using "./" instead ".*", empty channels will be "inf"

%--------------------organize the data
sg=[ch_Master(:,sg_Master_order) ch_Slave(:,sg_Slave_order)]; %ordered strain gages
ch_Master(:,sg_Master_order)=[];
ch_Master_names(sg_Master_order)=[];
ch_Slave(:,sg_Slave_order)=[];
ch_Slave_names(sg_Slave_order)=[];

ch_ex=[ch_Master,ch_Slave];
clear ch_Master;
clear ch_Slave;
ch_ex_names=[ch_Master_names;ch_Slave_names];

%--------Take out all 'None' channels
k=1;
for j=1:length(ch_ex_names)
    if strcmp(ch_ex_names{j},'None')
        index_None(k)=j;
        k=k+1;
    end
end
ch_ex(:,index_None)=[];
ch_ex_names(index_None)=[];

%------------sort the strains with ascending order of x_sg. x=0 is where the stoper is.
sort_index=exp_details.axis_sort_index;
U1=sg(:,1:3:end)*1000;%[mStrain]
U1=U1(:,sort_index);
U2=sg(:,2:3:end)*1000;%[mStrain]
U2=U2(:,sort_index);
U3=sg(:,3:3:end)*1000;%[mStrain]
U3=U3(:,sort_index);

% %-----correct shear sensitivity
%  gV=[0,0.1,0.95,-0.08];
%  [U1,U2,U3]=calc_shear_sensitivity4(U1,U2,U3,gV);
%  disp('shear sensitivity was corrected' )
% % % %------

[Sxx,Syy,Sxy,Uxx,Uyy,Uxy,note]=calculate_stress_strain_PlaneOption(U1,U2,U3,sg_angle,PlaneFlag);

acq132=struct('x_sg',x_sg,'y_sg',y_sg,'sg_angle',sg_angle,'U1',U1,'U2',U2,'U3',U3,'Uxx',Uxx,'Uyy',Uyy,'Uxy',Uxy,'Sxx',Sxx,'Syy',Syy,'Sxy',Sxy,'note',note);

% %-----The host of each sg
hostM(1:16)={'acq93_1'};
hostM(17:32)={'acq93_2'};
hostS(1:16)={'acq94_1'};
hostS(17:32)={'acq94_2'};
host=[hostM(sg_Master_order) hostS(sg_Slave_order)];
host=host(1:3:end);
host=host(sort_index);
acq132.host=host;

%--------Single strain gages (Not Rosettes)
%---find all channels with single strain gages
k_xx=1;
k_yy=1;
for j=1:length(ch_ex(1,:))
    if(strcmp(ch_ex_names{j},'Uxx'))
        Uxx_index(k_xx)=j;
        k_xx=k_xx+1;
    elseif(strcmp(ch_ex_names{j},'Uyy'))
        Uyy_index(k_yy)=j;
        k_yy=k_yy+1;
    end
end

if( exist('Uxx_index') )
    s.Uxx=ch_ex(:,Uxx_index)*1000;
    s.Uyy=ch_ex(:,Uyy_index)*1000;
    ch_ex(:,[Uxx_index Uyy_index])=[];
    ch_ex_names([Uxx_index Uyy_index])=[];
    s.x_Uxx=exp_details.ch_ex_x(Uxx_index);
    s.x_Uyy=exp_details.ch_ex_x(Uyy_index);
    acq132.s=s;
end

%---Extra channels

index=min([40 length(ch_ex(:,1))]);
for j=1:length(ch_ex(1,:))
    
    if(~ strcmp('None',ch_ex_names{j}))
        
        if( strfind(ch_ex_names{j},'ds'))
            if(strcmp(ch_ex_names{j}(1:2),'ds'))
                ch_ex(:,j)=ds_v2mu(ch_ex(:,j),ch_ex_names{j}(3:end));
                %ch_ex(:,j)=ch_ex(:,j)-mean(ch_ex(1:index,j)); % %Substruct the initial values
            elseif(strcmp(ch_ex_names{j}(1:3),'Ods'))
                filename=[ch_ex_names{j} '_cal.mat'];
                if exist([cal_dir '\' filename],'file')~=0
                    ch_ex(:,j)=Odsv_2mu([cal_dir '\' filename],ch_ex(:,j));
                    % ch_ex(:,j)=ch_ex(:,j)-mean(ch_ex(1:index,j)); % %Substruct the initial values
                end
            end
        end
        acq132=setfield(acq132,ch_ex_names{j},ch_ex(:,j));
    end
end





