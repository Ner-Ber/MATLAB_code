function DataStruct = phedi_getDataWithoutBigPic(exp_dir, eventNum, varargin)
    
    %% parameters
    defCell = {...
        'res',                  62500;...
        'lineNum',              'all';...
        'BlurROT',              1;...
        'UsrRes',               nan;...
        'numOOM',               3;...
        'slopePenalty',         0;...
        'expandLowBound',       -1;...
        'expandHighBound',      0;...
        'followMethod',         'RMS_minimization';...
        'phediRelative2first',  0;...
        'coarseShift',          0;...
        'preRowsTime',          3e-3;...
        'postRowsTime',         3e-3;...
        'prePhediTime',         1e-3;...
        'postPhediTime',        3e-4;...
        'smooth4Phedi',         0.001;...
        'PhediSafetyDist',      30;...
        'QuickAndDirty',        0;...
        'scratchWidth',         100*1e-6;...
        'scratchRegionMeters',  [0.07 0.09];...
        'defPhotoLoc',          0.08;...
        'Cr',                   1255;...
        'CohsvModl_h',          1e-7;...
        'skipPrecursor',        1;...
        'frontFallDetermn',     2;...
        'kind',                 'Slow';...
        'Nsmt',                 4;...
        'DoPlot',               0
        };
    Param = setDefaults4function_byName(defCell,varargin);
    %% displacement
    [UxStruct] = phedi_getUxWithoutBigPic(exp_dir, eventNum,Param);
    
    %% strains gages
    kind = Param.kind;
    exp_details = expDetailsRead(exp_dir);
    disp('Reading strain gages data...');
    sgDataStruct=acq132_event_get_data...
        (exp_dir,eventNum,'start','end',1,'Uxx','Uyy','Uxy','Sxx','Syy','Sxy','U1','U2','U3','x_sg','y_sg','sg_angle','F','N');
    %---Fix for shear sensitivity (this part is taken from 'anls_front_in_space.m')
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
    
    
    %% phedi velocity
    TemporalVecPix_D = diff(UxStruct.TemporalVecPix);
    TemporalVecSec_D = diff(UxStruct.TemporalVecSec);
    Vx_PixPrFrm_preSmth = bsxfun(@rdivide,diff(UxStruct.PhediLocationPixShifted,1,1),TemporalVecPix_D(:));
    Vx_MetPrSec_preSmth = bsxfun(@rdivide,diff(UxStruct.PhediLocationPixShifted,1,1)/UxStruct.res,TemporalVecSec_D(:));
    
    
    %% smooth and get velocity
    Nsmt = Param.Nsmt;
    PhediLocation = conv2(UxStruct.PhediLocationPixShifted/UxStruct.res,ones(Nsmt,1)/Nsmt,'same');
    PhediLocation([1:Nsmt,end-Nsmt:end],:) = nan;
    
    PhediLocationPix = conv2(UxStruct.PhediLocationPixShifted,ones(Nsmt,1)/Nsmt,'same');
    PhediLocationPix([1:Nsmt,end-Nsmt:end],:) = nan;
    
    Vx_PixPrFrm = bsxfun(@rdivide,diff(PhediLocationPix,1,1),TemporalVecPix_D(:));
    Vx_MetPrSec = bsxfun(@rdivide,diff(PhediLocation,1,1),TemporalVecSec_D(:));
    
    time4vel_sec = movmean(UxStruct.TemporalVecSec,2,'endpoints','discard');
    time4vel_frm = movmean(UxStruct.TemporalVecPix,2,'endpoints','discard');
    
    meanPhedLocM = mean(PhediLocation,2);
    meanPhedVelMpS = mean(Vx_MetPrSec,2);
    
    %% save stuff
    DataStruct.UxStruct = UxStruct;
    DataStruct.UxStruct.Vx_PixPrFrm_preSmth = Vx_PixPrFrm_preSmth;
    DataStruct.UxStruct.Vx_MetPrSec_preSmth = Vx_MetPrSec_preSmth;
    DataStruct.UxStruct.PhediLocation = PhediLocation;
    DataStruct.UxStruct.PhediLocationPix = PhediLocationPix;
    DataStruct.UxStruct.Vx_PixPrFrm = Vx_PixPrFrm;
    DataStruct.UxStruct.Vx_MetPrSec = Vx_MetPrSec;
    DataStruct.UxStruct.time4vel_sec = time4vel_sec;
    DataStruct.UxStruct.time4vel_frm = time4vel_frm;
    DataStruct.UxStruct.meanPhedVelMpS = meanPhedVelMpS;
    DataStruct.UxStruct.meanPhedLocM = meanPhedLocM;
    
    
    DataStruct.sgDataStruct = sgDataStruct;
    
    
    
    
end