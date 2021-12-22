function acq132=acq132_convert_raw_to_data1(exp_dir,ch_093,ch_094)

% uses "exp_details.txt" to convert units , correct offset and s.g rotation, calculate stress and strain.
%It makes the final data
%exp_dir is needed to use event1 for offset substruction

%------------Read exp_details.txt
exp_details=expDetailsRead(exp_dir);
factor_093=exp_details.factor_093;
offset_093=exp_details.offset_093;
offset_094=exp_details.offset_094;
factor_094=exp_details.factor_094;
sg_093_order=exp_details.sg_093_order;
sg_094_order=exp_details.sg_094_order;
sg_angle=exp_details.sg_angle;
x_sg=exp_details.x_sg;%notice the positive x axis on the exp_details is in the pushing direction.Here we work in the camera frame (opposite direction) thus must be sorted
ch_093_names=exp_details.ch_093_names;
ch_094_names=exp_details.ch_094_names;
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
    path_093=[cal_dir '\acq132_093\multivent'];
    path_094=[cal_dir '\acq132_094\multivent'];
else
    path_093=[exp_dir '\acq132_093\multivent'];
    path_094=[exp_dir '\acq132_094\multivent'];
end

ch_offset=acq132_event_read(path_093,1);%using event1 for offset correction
ch_offset=mean(ch_offset);
ch_number_for_offset = find(offset_093==0);
ch_093(:,ch_number_for_offset)=ch_093(:,ch_number_for_offset)-repmat(ch_offset(ch_number_for_offset),length(ch_093(:,1)),1);%offset correction// comment this line for disable
ch_number_for_offset = find(offset_093~=0);
ch_093(:,ch_number_for_offset)=ch_093(:,ch_number_for_offset)-repmat(offset_093(ch_number_for_offset),length(ch_093(:,1)),1);%offset correction with specified values// comment this line for disable
ch_093=ch_093./repmat(factor_093,length(ch_093(:,1)),1);%unit conversion notice the using "./" instead ".*", empty channels will be "inf"

%--offset for ch_094_ex
ch_offset=acq132_event_read(path_094,1);%using event1 for offset correction
ch_offset=mean(ch_offset);
ch_number_for_offset = find(offset_094==0);
ch_094(:,ch_number_for_offset)=ch_094(:,ch_number_for_offset)-repmat(ch_offset(ch_number_for_offset),length(ch_094(:,1)),1);%offset correction// comment this line for disable
ch_number_for_offset = find(offset_094~=0);
ch_094(:,ch_number_for_offset)=ch_094(:,ch_number_for_offset)-repmat(offset_094(ch_number_for_offset),length(ch_094(:,1)),1);%offset correction with specified values// comment this line for disable
ch_094=ch_094./repmat(factor_094,length(ch_094(:,1)),1);%unit conversion notice the using "./" instead ".*", empty channels will be "inf"

%--------------------organize the data
sg=[ch_093(:,sg_093_order) ch_094(:,sg_094_order)]/2; %ordered strain gages, factor 1/2 is needed for sg's gage factor
N=ch_093(:,16);%N=ch_094(:,16);
F=ch_093(:,32);
[~,column]=find(factor_094(1,16:end)~=0);%All extra channels ,not sg and loads.All data with factor=0 is ignored.
ch_094_ex=ch_094(:,15+column);
clear ch_094 ch_093; %the raw data isn't used later

% ---------convert units and cal to ch_094_ex ------
%---not automatic choose channels by hand (in future maybe useful to make
%mfile per experimemt saved at the exp dirictory that treats all external
%channels in different way

if ~isempty(column)
    %     ds_offset=2.5;%analog voltage offset
    % ch_094_ex(:,3)=ch_094_ex(:,3)+ds_offset;
    %     % ds_offset=1;%analog voltage offset
    % ch_094_ex(:,[2 4:end])=ch_094_ex(:,[2 4:end])+ds_offset;
    
    
    %------ Ods & ds calibration
    index=min([40 length(ch_094_ex(:,1))]);
    
    if length(ch_094_ex(1,:))>1
        
        filename='Odscal_ch_094_ex2.mat';
        if exist([cal_dir '\' filename],'file')~=0
            ch_094_ex(:,2)=Odsv_2mu([cal_dir '\' filename],ch_094_ex(:,2));
            ch_094_ex(:,2)=ch_094_ex(:,2)-mean(ch_094_ex(1:index,2)); % %Substruct the initial values
        else
%              ch_094_ex(:,2)=ds_v2mu(ch_094_ex(:,2),2507);
%              ch_094_ex(:,2)=ch_094_ex(:,2)-mean(ch_094_ex(1:index,2)); % %Substruct the initial values
        end
        
        filename='Odscal_ch_094_ex3.mat';
        if exist([cal_dir '\' filename],'file')~=0 %correction for offset was done
            ch_094_ex(:,3)=Odsv_2mu([cal_dir '\' filename],ch_094_ex(:,3));
            ch_094_ex(:,3)=ch_094_ex(:,3)-mean(ch_094_ex(1:index,3)); % %Substruct the initial values
            
        else
            %             ch_094_ex(:,3)=ds_v2mu_2508(ch_094_ex(:,3));
            %             ch_094_ex(:,3)=ch_094_ex(:,3)-mean(ch_094_ex(1:index,3)); % %Substruct the initial values
            %             ch_094_ex(:,3)=-ch_094_ex(:,3);
        end
        
        filename='Odscal_ch_094_ex4.mat';
        if exist([exp_dir '\cal\' filename],'file')~=0
            ch_094_ex(:,4)=Odsv_2mu([cal_dir '\' filename],ch_094_ex(:,4));
            ch_094_ex(:,4)=ch_094_ex(:,4)-mean(ch_094_ex(1:index,4)); %Substruct the initial values
            
        else
            %ch_094_ex(:,4)=ds_v2mu_2404(ch_094_ex(:,4));
            %ch_094_ex(:,4)=ch_094_ex(:,4)-mean(ch_094_ex(1:index,4)); %Substruct the initial values
            ch_094_ex(:,4)=-ds_v2mu(ch_094_ex(:,4),2508);
            ch_094_ex(:,4)=ch_094_ex(:,4)-mean(ch_094_ex(1:index,4)); % %Substruct the initial values

        end
        
        filename='Odscal_ch_094_ex5.mat';
        if exist([cal_dir '\' filename],'file')~=0
            ch_094_ex(:,5)=Odsv_2mu([cal_dir '\' filename],ch_094_ex(:,5));
            ch_094_ex(:,5)=ch_094_ex(:,5)-mean(ch_094_ex(1:index,5)); % %Substruct the initial values
        else
            %sprintf('no calibration file for ch_094_ex5')
        end
        
        
        filename='Odscal_ch_094_ex6.mat';
        if exist([exp_dir '\cal\' filename],'file')~=0
            ch_094_ex(:,6)=Odsv_2mu(exp_dir,filename,ch_094_ex(:,6));
            ch_094_ex(:,6)=ch_094_ex(:,6)-mean(ch_094_ex(1:index,6)); % %Substruct the initial values
        else
            %         ch_094_ex(:,6)=ds_v2mu_1990fs(ch_094_ex(:,6));
            %         ch_094_ex(:,6)=ch_094_ex(:,6)-mean(ch_094_ex(1:index,6)); % %Substruct the initial values
        end
        
    end
end

%--------------------calculate strain stress
%-in the Camera axis frame (positive is oposite to pushing direction)

Uxy=-(sg(:,1:3:end)-sg(:,3:3:end))/2;  %calculated shear strain %the - is needed to change the direction of x axis. previosly in the pushing directio, now opposite
Uyy=sg(:,2:3:end);                     %calculated normal strain
Uxx=(sg(:,1:3:end)+sg(:,3:3:end))-Uyy; %calculated transvers strain

%------------sort the strains with ascending order of x_sg. x=0 is where the stoper is.
sort_index=exp_details.axis_sort_index;

Uxy=Uxy(:,sort_index)*1000; %[mStrain]
Uyy=Uyy(:,sort_index)*1000; %[mStrain]
Uxx=Uxx(:,sort_index)*1000; %[mStrain]

%------------------- rotate the strain (correction)
for j=1:length(Uxx(1,:))
    u(1,1,:)=Uxx(:,j);
    u(1,2,:)=Uxy(:,j);
    u(2,1,:)=Uxy(:,j);
    u(2,2,:)=Uyy(:,j);
    R=[cos(sg_angle(j)) sin(sg_angle(j)); -sin(sg_angle(j)) cos(sg_angle(j))]; % rotation => u=R(theta)u'R(-theta)
    u=multiprod(u,R);
    R=[cos(sg_angle(j)) -sin(sg_angle(j)); sin(sg_angle(j)) cos(sg_angle(j))];
    u=multiprod(R,u);
    Uxx(:,j)=u(1,1,:);
    Uxy(:,j)=u(1,2,:);
    Uyy(:,j)=u(2,2,:);
end

%---------now calc stresses using hook law and assuming Szz=0
sigma=1/3;%poisson ratio
E=3E+9/10^6 ;%youngs modulus [MPa]
E=E/1000; %because U is in mStrain

Sxx=-E/(1-sigma^2)*(Uxx+sigma*Uyy);  % "-" is added positive stress is compresion
Syy=-E/(1-sigma^2)*(Uyy+sigma*Uxx); % "-" is added positive stress is compresion
Sxy=E/(1+sigma)*Uxy;


acq132=struct('x_sg',x_sg,'Uxx',Uxx,'Uyy',Uyy,'Uxy',Uxy,'Sxx',Sxx,'Syy',Syy,'Sxy',Sxy,'N',N,'F',F,'ch_094_ex',ch_094_ex);

%Uzz=-sigma/(1-sigma)*(Uxx+Uyy);
%Sxx=E/((1+sigma)*(1-2*sigma))*((1-sigma)*Uxx+sigma*(Uyy+Uzz));
%Syy=E/((1+sigma)*(1-2*sigma))*((1-sigma)*Uyy+sigma*(Uxx+Uzz));
%Sxy=E/(1+sigma)*Uxy;




