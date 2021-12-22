
function [D1,D2,T,V1,V2,V3,V_rupture] = Find_movement_of_sensors(new_calibration_file,event_file,fileType)
  
D_sensors=0.045;    %%the measured distance between the sensors
    %% get data from the new calibration file
%     A=xlsread(new_calibration_file);

    v1 = new_calibration_file(1:end,1);
    v2 = new_calibration_file(1:end,2);
    d = new_calibration_file(1:end,3);
    
    %% get data from the event file
    if strcmp(fileType, 'csv')
        start=19;
        B=xlsread(event_file);
        T=B(start:end,1);        %%time 
        T=T-T(1);
        V1=B(start:end,2);    %%left sensor
        V2=B(start:end,3);    %%right sensor
        V3=B(start:end,4);    %%calibration sensor
    elseif strcmp(fileType, 'mat')
        event_structure = load(event_file);
        T = event_structure.Trace_1(:,1);
        T=T-T(1);
        V1 = event_structure.Trace_1(:,2);
        V2 = event_structure.Trace_2(:,2);
        V3 = event_structure.Trace_3(:,2);
    end
    
    %% find rupture velocity
    
    %# Defines the movement starting-point to be the first time in which...
    ...the measured value exceed the max value before the rupture
        
    max_V1=max(V1(1:round(0.4*length(V1))));
    max_V2=max(V2(1:round(0.4*length(V2))));
    t1=T(find(V1>max_V1,1));  
    t2=T(find(V2>max_V2,1));
    V_rupture=D_sensors/(t2-t1);

    %% find linear accent of calibration
    % the rational for this step is cutting off the "bad" measurments,
    % those who stick up into the area where the the calibrqation is not
    % deterministic. In order to do this the we need to first find the
    % right range, in order to find the right range the noise needs to be
    % smoothended out and the derivaive should be found. 

    %% use the new calibration to interpret the vectors D1,D2 of the distance at all time
    %--- create logicals for picking range for interpolation. this must be
    %done due to influence of  boundary consitions on the csaps function.
    rangeLogical1 = min(V1)<v1 & v1<max(V1);
    rangeLogical2 = min(V2)<v2 & v2<max(V2);
    
    %--- use logicals in order to interpolate only with the relevant part
    %of vi and d.
    D1=csaps(v1(rangeLogical1),d(rangeLogical1),0.3,V1);
    D2=csaps(v2(rangeLogical2),d(rangeLogical2),3,V2);
   
    %% plots
%     
%     %---plot1 - the measurment: volts vs time
%     figure
%     subplot(2,1,1)
%     hold all
%     plot(T,V1,'.','color','b')
%     plot(T,V2,'.','color','g')
%     plot(T,V3,'.','color','r')
%     hold off
%     ylim([-0.5 5])
%     title('\fontsize{16}Event Measurment');
%     xlabel('Time [s]');
%     ylabel('Voltage [v]');
%     legend('V1','V2','V3','Location','NorthEastOutside');
% 
%     %---plot2 - the movement: distance vs time
%     subplot(2,1,2)
    
    figure;
    
    hold all
    plot(T,D1,'.','color','b')
    plot(T,D2,'.','color','g')
    title({'\fontsize{16}Event Movement' ; ['\fontsize{12}V rupture=' num2str(V_rupture) '[m/s]']});
    xlabel('Time [s]');
    ylabel('Displacement [mm]');
    legend('V1','V2','V3','Location','NorthEastOutside');
    hold off
    
    

    
end