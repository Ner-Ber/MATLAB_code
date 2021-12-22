function [AvgPhediStruct_1_2_3, AvgPhediStruct_6_7_8, AvgPhediStruct_all] = plot_z_diff_comparison(z_diff_struct,PhediStructCell, Ev, Cf)
    
    figure; hold on;
    AvgPhediStruct_1_2_3 = phedi_averagePhedi(z_diff_struct.PhediStructMultipLines.line_1_2_3);
    Name='lines 1:3';
    x_Phedi_1_2_3 = AvgPhediStruct_1_2_3.X_mean;
    vel_phedi_1_2_3 = AvgPhediStruct_1_2_3.Vel_mean_x;
    plot(x_Phedi_1_2_3,vel_phedi_1_2_3,'DisplayName', Name)
    
    
    AvgPhediStruct_6_7_8 = phedi_averagePhedi(z_diff_struct.PhediStructMultipLines.line_6_7_8);
    Name='lines 6:8';
    x_Phedi_6_7_8 = AvgPhediStruct_6_7_8.X_mean;
    vel_phedi_6_7_8 = AvgPhediStruct_6_7_8.Vel_mean_x;
    plot(x_Phedi_6_7_8,vel_phedi_6_7_8,'DisplayName', Name)
    
    
    Name = 'all lines avg';
    AvgPhediStruct_all = phedi_averagePhedi(PhediStructCell{Ev});
    x_Phedi = AvgPhediStruct_all.X_mean;
    vel_phedi = AvgPhediStruct_all.Vel_mean_x;
    plot(x_Phedi,vel_phedi,'DisplayName', Name)
    
    title(sprintf('Cf=%d, Ev=%d', round(Cf), Ev))
    
end