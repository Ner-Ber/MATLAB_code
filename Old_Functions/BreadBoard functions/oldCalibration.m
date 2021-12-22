%%%calibration_file='C:\Users\owner\Documents\Shahar_Neri\friction_videos\kiul_check\C6a.csv';
%%%half=4240;
function [V1,V2,V3,D,T] = oldCalibration(calibration_file,plot1,plot2,half)
    %% getting data from the early calibration of the "calibration sensor"...
    %taking only points from the peak and on to prevent double values for the same v
    filename='C:\Users\owner\Documents\Shahar_Neri\displacement sensors calibration.xlsx';
    A=xlsread(filename);
    peak=87;   
    d=A(peak:end,13);
    v=A(peak:end,14);

    %% getting data from the calibration measurment
    B=xlsread(calibration_file);
    start=19;
    T=B(19:end,1);   %%time 
    T=T-T(1);
    V1=B(start:end,2);    %%left sensor
    V2=B(start:end,3);    %%right sensor
    V3=B(start:end,4);    %%calibration sensor
    half=half-start;
    
    %% use the early calibration to interpret the vector D of the distance at all time
    D=zeros(size(V3));
    for i=1:size(V3)
        D(i)=interp1(v,d,V3(i));
    end

    %% plot
    switch plot1
        case 1
            figure
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
    end
    
    switch plot2
        case 1
            figure
            hold all
            plot(D,V1,'.','color','b')
            plot(D,V2,'.','color','g')
            hold off
            ylim([-0.5 5])
            title('\fontsize{16}Calibration');
            xlabel('Displacement [mm]');
            ylabel('Voltage [v]');
            legend('V1','V2','Location','NorthEastOutside');
        case 2
            figure
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
            legend('V1','V1a','V2','V2a','Location','NorthEastOutside');
            
    end
end