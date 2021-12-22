% cutSignalsFromReturningWaves will cut the signals which will be presented when a
% returning wave meets the front.
% the returning wave-speeed is determined by usr. 



function [returnLogical, returnIdx_t] = cutSignalsFromReturningWaves(WaveSpeed,PhediStruct, AvgPhediStruct )
    
    %--- calc time for some wave to get to from the scratches to the boundary and back
    WaveTrvlDist = 2*(0.15-PhediStruct.PhediData.PhotoLocation);
    WaveTrvlTime = WaveTrvlDist/WaveSpeed;

    %--- maybe with a reflected front from the boundary? (doesn't seem to work well)
%     [~,I] = min(abs(PhediStruct.BigPicRotStruct.frontVelLoc_interpM - 0.14));
%     WaveSpeed = PhediStruct.BigPicRotStruct.frontVel_interpMperS(I);
%     WaveTrvlDist = abs(0.15-PhediStruct.PhediData.PhotoLocation)+abs(0.14-0.15);
%     WaveTrvlTime = WaveTrvlDist/WaveSpeed;
    
    x_Phedi = AvgPhediStruct.X_mean;
    vel_phedi = AvgPhediStruct.Vel_mean_x;
    
    [~,returnIdx_t] = min(abs(AvgPhediStruct.T_mean-WaveTrvlTime));
    
    [waveReturn_x, ~] = signal_ChangeTime2Space_mapping(WaveTrvlTime, PhediStruct.PhediData.PhotoLocation, PhediStruct.BigPicRotStruct);
    
    returnLogical = x_Phedi>=waveReturn_x;
    
    
%     returnLogical = x_Phedi>=AvgPhediStruct.X_mean(returnIdx_t);
    
    
    
end