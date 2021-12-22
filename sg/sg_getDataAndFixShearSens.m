function DataStruct_wSG = sg_getDataAndFixShearSens(DataStruct,varargin)
% sgDataStruct = sg_getDataAndFixShearSens(DataStruct,varargin)
%
%sg_getDataAndFixShearSens is intended to work in the function
%phedi_createStructureWAdd and add data from the strain gages (shear
%sensitivity fixed) to the main DataStruct.
%
% VARARGIN: (DEFAULT)
% kind: front type, can be 'Slow','Rayleigh','Super' (def: 'Super')

%% set defaults and get fields
kind = setDefaults4function(varargin,'Slow');
BigPicRotStruct = DataStruct.BigPicRotStruct;
eventNum = DataStruct.ExperimentData.eventNum;
exp_dir = DataStruct.ExperimentData.ExpHour;
exp_details = DataStruct.ExperimentData;
%% get sg data
disp('Reading strain gages data...');
sgDataStruct=acq132_event_get_data...
    (exp_dir,eventNum,'start','end',1,'Uxx','Uyy','Uxy','Sxx','Syy','Sxy','U1','U2','U3','x_sg','y_sg','sg_angle','F','N','trigger','note','host');
%--- add velocities and vectors for each sg and update DataStruct:
sgDataStruct = sg_createVectorsForSg(sgDataStruct,BigPicRotStruct);


%% Fix for shear sensitivity (this part is taken from 'anls_front_in_space.m')
% acqE=acq132_event_get_data(exper,event,'start','end',smtAcq,'U1','U2','U3','x_sg','y_sg','sg_angle','host');
%--- Shift U1&U3
sub_sample_index=100;
dt=(sgDataStruct.t(2)-sgDataStruct.t(1))/sub_sample_index;%(msec)
if(strcmp(kind,'Super'))
    shift_t=1.6e-4;%0.3mm/2200m/s=0.14mus
else
    shift_t=2.5e-4;%0.3mm/1200m/s=0.25mus
end

d=ceil(shift_t/dt);

t_spline=(sgDataStruct.t(1):dt:sgDataStruct.t(end));
for j=1:length(sgDataStruct.x_sg)
    if (exp_details.sg_angle(j)==0) %Otherwise some thing more elaborated should be done.
        U1_spline=spline(sgDataStruct.t,sgDataStruct.U1(:,j),t_spline)';
        U2_spline=spline(sgDataStruct.t,sgDataStruct.U2(:,j),t_spline)';
        U3_spline=spline(sgDataStruct.t,sgDataStruct.U3(:,j),t_spline)';
        
        U1_spline=circshift(U1_spline,d);
        U3_spline=circshift(U3_spline,-d);
        
        sgDataStruct.U1(:,j)=U1_spline(1:sub_sample_index:end);
        sgDataStruct.U2(:,j)=U2_spline(1:sub_sample_index:end);
        sgDataStruct.U3(:,j)=U3_spline(1:sub_sample_index:end);
    end
end
disp(['shift=' num2str(shift_t)] )

sgDataStruct.gV=[0,0.1,0.95,-0.08];
%sgDataStruct.gV=[0.0,0.15,0.95,-0.08];
%sgDataStruct.gV=[0,0,1,0];
[sgDataStruct.U1,sgDataStruct.U2,sgDataStruct.U3]=calc_shear_sensitivity4(sgDataStruct.U1,sgDataStruct.U2,sgDataStruct.U3,sgDataStruct.gV);
[sgDataStruct.Sxx,sgDataStruct.Syy,sgDataStruct.Sxy,sgDataStruct.Uxx,sgDataStruct.Uyy,sgDataStruct.Uxy,~]=calculate_stress_strain_PlaneOption(sgDataStruct.U1,sgDataStruct.U2,sgDataStruct.U3,sgDataStruct.sg_angle);
disp('shear sensitivity was corrected' )

%% add to DataStruct
DataStruct_wSG = DataStruct;
DataStruct_wSG.SgData = sgDataStruct;
end