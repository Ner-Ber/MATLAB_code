function acq132=mixVK_acq132_convert_raw_to_data(exp_dir,ch_Master,ch_Slave)

% uses "exp_details.txt" to convert units , correct offset and s.g rotation, calculate stress and strain.
%It makes the final data
%exp_dir is needed to use event1 for offset substruction
%specific routine for the block 200mm with Vishay and Kulite. 
%Master = Kulite, Slave=Vishay.

%------------Read exp_details.txt
exp_details=mixVK_expDetailsRead(exp_dir);
factor_Master=exp_details.factor_Master;
offset_Master=exp_details.offset_Master;
offset_Slave=exp_details.offset_Slave;
factor_Slave=exp_details.factor_Slave;
sg_Master_order=exp_details.sg_Master_order;
sg_Slave_order=exp_details.sg_Slave_order;
x_sg=exp_details.x_sg;
ch_Master_names=exp_details.ch_Master_names;
ch_Slave_names=exp_details.ch_Slave_names;
MasterVref=exp_details.MasterVref;
MasterCurrent=exp_details.MasterCurrent*1e-6;
MasterResistances=exp_details.MasterResistances;
%Calculation of gains of Master (Kulite)
factor_Master([1:15 17:31])=factor_Master([1:15 17:31]).*MasterCurrent;% it has to be multiplied by the gain (NL) and resistances.
%factor_Master([1:15 17:31])=factor_Master([1:15 17:31]).*MasterResistances*MasterCurrent*1e-6;% it has to be multiplied by the gain (NL).
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
% Calculation of resistances from the offset values
%MasterResistances(1:16)=(ch_offset(1:16)-MasterVref(1))./factor_Master(1:16); MasterResistances(17:32)=(ch_offset(17:32)-MasterVref(2))./factor_Master(17:32);
%MasterResistances=[1304 1323 1275 1393 1373 1434 1284 1265 1245 1180 1156 1171 1206 1203 1203 0 1394 1328 1347 1386 1312 1354 1184 1172 1171 1108 1158 1136 1203 1128 1177 0];
%factor_Master([1:15 17:31])=factor_Master([1:15 17:31]).*MasterResistances([1:15 17:31]);% it has to be multiplied by the gain (NL) and resistances.
factor_Master([1:15 17:31])=factor_Master([1:15 17:31]).*MasterResistances;

ch_number_for_offset = find(offset_Master==0);
ch_Master(:,ch_number_for_offset)=ch_Master(:,ch_number_for_offset)-repmat(ch_offset(ch_number_for_offset),length(ch_Master(:,1)),1);%offset correction// comment this line for disable
ch_number_for_offset = find(offset_Master~=0);
ch_Master(:,ch_number_for_offset)=ch_Master(:,ch_number_for_offset)-repmat(offset_Master(ch_number_for_offset),length(ch_Master(:,1)),1);%offset correction with specified values// comment this line for disable
ch_Master=ch_Master./repmat(factor_Master,length(ch_Master(:,1)),1);%unit conversion notice the using "./" instead ".*", empty channels will be "inf"

%Calculation of the gain 
 %g=nlGain(ch_Master(:,[1:15 17:31]));
 gain=132*ones(length(ch_Master(:,1)),32);
 %gain(:,[1:15 17:31])=g;
% gain(:,[1:15 17:31])=[105.7665   74.5891   96.3485  112.7763  105.3140  111.3248  106.1466  103.2258  107.4539  118.9276  106.2668  105.0868  121.8213   98.9873   96.7735  132.2041  122.7879  106.9132  140.2981 103.4137  109.0602  120.0223  117.1822  115.3758  137.5575  103.5922   94.2662  161.5935  113.8072   93.9353];

gain(:,[16 32])=1;
%07/08/2014
%  ch_Master=ch_Master./repmat(gain,length(ch_Master(:,1)),1);
 ch_Master=ch_Master./gain;

% ch_Master=ch_Master./repmat(factor_Master,length(ch_Master(:,1)),1);%unit conversion notice the using "./" instead ".*", empty channels will be "inf"

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

[Sxx,Syy,Sxy,Uxx,Uyy,Uxy,note]=ang_mixVK_calculate_stress_strain(sg,exp_details);
acq132=struct('x_sg',x_sg,'Uxx',Uxx,'Uyy',Uyy,'Uxy',Uxy,'Sxx',Sxx,'Syy',Syy,'Sxy',Sxy,'note',note,'sg',sg);

%---Extra channels

index=min([40 length(ch_ex(:,1))]);
for j=1:length(ch_ex(1,:))
    
    if(~ strcmp('None',ch_ex_names{j}))
        
        if( strfind(ch_ex_names{j},'ds'))
            if(strcmp(ch_ex_names{j}(1:2),'ds'))
                ch_ex(:,j)=ds_v2mu(ch_ex(:,j),ch_ex_names{j}(3:end));
                ch_ex(:,j)=ch_ex(:,j)-mean(ch_ex(1:index,j)); % %Substruct the initial values
            elseif(strcmp(ch_ex_names{j}(1:3),'Ods'))
                filename=[ch_ex_names{j} '_cal.mat'];
                if exist([cal_dir '\' filename],'file')~=0
                    ch_ex(:,j)=Odsv_2mu([cal_dir '\' filename],ch_ex(:,j));
                    ch_ex(:,j)=ch_ex(:,j)-mean(ch_ex(1:index,j)); % %Substruct the initial values
                end
            end
        end
        acq132=setfield(acq132,ch_ex_names{j},ch_ex(:,j));
    end
end





