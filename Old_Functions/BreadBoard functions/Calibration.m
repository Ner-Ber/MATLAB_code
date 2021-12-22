
% calibration_file='C:\Users\owner\Documents\Shahar_Neri\friction_videos\19.9\12\C12p.csv';

function [V1,V2,V3,D,T,new_calibration_file] = Calibration(calibration_file,fileType)
%fileType - should be written either 'mat', or 'csv' to fit the calibration
%file type enterred.

%% choose computer
computer='\\ELSA\c'; %% Elsa
%computer='C:';  %% This_computer
    
    %%
    [filesFolder, calibration_file_name, ~] = fileparts(calibration_file);

    % [calibration_file,half]=Choose_calibration(calibration_number);
        %% getting data from the early calibration of the "calibration sensor"...
        %taking only points from the peak and on to prevent double values for the same v
        filename=[computer '\Users\owner\Documents\Shahar_Neri\displacement sensors calibration.xlsx'];
        A=xlsread(filename);
        peak=87;   
        d=A(peak:end,13);
        v=A(peak:end,14);
        
        %% getting data from the calibration measurment
        if strcmp(fileType, 'csv')
            B=xlsread(calibration_file);
            start=19; %%old scope
            s=1;  %%old scope
    %         start=23;  %%new scope
    %         s=2;       %%new scope
            T=B(start:end,s);   %%time 
            T=T-T(1);
            V1=B(start:end,s+1);    %%left sensor
            V2=B(start:end,s+2);    %%right sensor
            V3=B(start:end,s+3);    %%calibration sensor
            
        elseif strcmp(fileType, 'mat')
            
            calibration_structure = load(calibration_file);
            T = calibration_structure.Trace_1(:,1);
            V1 = calibration_structure.Trace_1(:,2);
            V2 = calibration_structure.Trace_2(:,2);
            V3 = calibration_structure.Trace_3(:,2);
            start = 0;
            
        end

        %% finding the middle of the calibration to avoid double values
        half=find(V3==max(V3));
        half=half(1);        %% in case there are several places with max value
        half=half-start;
        
        %% use the early calibration to interpret the vector D of the distance at all time
        D=interp1(v,d,V3);
        
        %% plot
    %     switch plot1
    %         case 1
                figure
                subplot(2,1,1)
                hold all
                plot(T,V1,'.','color','b')
                plot(T,V2,'.','color','g')
                plot(T,V3,'.','color','r')
                hold off
                ylim([-0.5 5])
                title('\fontsize{16}Measurment');
                xlabel('Time [s]');
                ylabel('Voltage [v]');
                legend('V1','V2','V3','Location','NorthEastOutside');
    %     end

    %     switch plot2
    %         case 1
    %             figure
    %             hold all
    %             plot(D,V1,'.','color','b')
    %             plot(D,V2,'.','color','g')
    %             hold off
    %             ylim([-0.5 5])
    %             title('\fontsize{16}Calibration');
    %             xlabel('Displacement [mm]');
    %             ylabel('Voltage [v]');
    %             legend('V1','V2','Location','NorthEastOutside');
    %         case 2
                subplot(2,1,2)
                hold all
                plot(D(1:half),V1(1:half),'.','color',[0 0 1])
                plot(D(half:end),V1(half:end),'.','color',[0.4 0.4 1])
                plot(D(1:half),V2(1:half),'.','color',[0 1 0])
                plot(D(half:end),V2(half:end),'.','color',[0.4 1 0.4])
                hold off
                ylim([-0.5 5])
                title('\fontsize{16}Calibration');
                xlabel('Displacement [mm]');
                ylabel('Voltage [v]');
                legend('V1a','V1b','V2a','V2b','Location','NorthEastOutside');

    %     end

    %% smooth
    % Calib_V1=smooth(V1a,10);
    % Calib_V2=smooth(V2a,10);
    % Displacement=smooth(Da,30);

    %% save the new calibration for V1 & V2

    % new_calibration_file=[calibration_file '_new.xlsx'];
    V1a=V1(1:half);
    V2a=V2(1:half);
    Da=D(1:half);
    new_calibration_file = [V1a,V2a,Da];
    
    
    save(fullfile(filesFolder, [calibration_file_name,'_new.mat']),...
        'new_calibration_file');
  


end
