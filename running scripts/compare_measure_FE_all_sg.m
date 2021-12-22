
%sg_getDataAndFixShearSens is intended to work in the function
%phedi_createStructureWAdd and add data from the strain gages (shear
%sensitivity fixed) to the main DataStruct.
%
% VARARGIN: (DEFAULT)
% kind: front type, can be 'Slow','Rayleigh','Super' (def: 'Super')

%% general
exper = my_dir;
expArray = struct;
for exp_idx=1:12
    exp_dir = exper{exp_idx};
    outdata = phantomReadMeta([exp_dir,'\Ph']);
    GammaMat = nan(30,13);
    CfMat = nan(outdata.NumEvents,15);
    for evIdx = 1:outdata.NumEvents
        try
            eventNum = evIdx;
            sg_vec = 1:15;
            %% big pic struct
            BigPicRotStruct = phantom_BigPicStruct(exp_dir, eventNum,5e-3,1e-10,5e-3,1,1,'all',2);
            
            if ~phantom_isPrecursor(BigPicRotStruct)
                %% set defaults and get fields
                kind = 'Slow';
                exp_details = expDetailsRead(exp_dir);
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
                
                
                %% update structure
                updated_sgDataStruct = sg_createVectorsForSg(sgDataStruct,BigPicRotStruct);
                
                
                %% space vectors for sg
                CfSmooth = 10;
                x_mins_x_tip_sg = zeros(size(updated_sgDataStruct.t_mins_t_tips));
                for i = 1:size(updated_sgDataStruct.t_mins_t_tips,2)
                    try
                        x_min_Ttip_sgi = signal_ChangeTime2Space_mapping(...
                            updated_sgDataStruct.t_mins_t_tips(:,i), updated_sgDataStruct.x_sg(i)*1e-3, BigPicRotStruct);
                    catch
                        disp();
                    end
                    x_mins_x_tip_sg(:,i) = x_min_Ttip_sgi(:);
                end
                %--- backup old and save new
                updated_sgDataStruct.x_mins_x_tips_constCf = updated_sgDataStruct.x_mins_x_tips;
                updated_sgDataStruct.x_mins_x_tips = x_mins_x_tip_sg;
                
                
                for s_i = sg_vec
                    %% find fracture energy
                    solAtSG = gamma_findGammaFit_ByAmpUxx(updated_sgDataStruct,BigPicRotStruct, s_i);
                    GammaMat(evIdx,s_i) = solAtSG.Gamma;
                    CfMat(evIdx,s_i) = solAtSG.v;
                end
                
                %% plot
                % Colors = MyVaryColor(length(sg_vec));
                % figure; hold on;
                % for s_i = sg_vec
                %     %% find fracture energy
                %     solAtSG = gamma_findGammaFit_ByAmpUxx(updated_sgDataStruct,BigPicRotStruct, s_i);
                %     %%
                %     %--- plot Uxx
                %     plot(updated_sgDataStruct.x_mins_x_tips(:,s_i),1e-3*(updated_sgDataStruct.Uxx(:,s_i)-mean(updated_sgDataStruct.Uxx(1:100,s_i))),...
                %         '.-','Color',Colors(sg_vec==s_i,:));
                %     plot(solAtSG.x,solAtSG.Uxx,...
                %         'Color',Colors(sg_vec==s_i,:));
                % end
            end
        catch
            display(['ev ',num2str(evIdx),' in exp ',exper{exp_idx},' failed']);
        end
    end
    expArray(exp_idx).CfMat = CfMat;
    expArray(exp_idx).GammaMat = GammaMat;
end