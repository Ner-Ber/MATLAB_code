%%%gets a serial number for a calibration file and gives his location and
%%%the middle point of the calibration (if its a back-and-forth
%%%calibration)
function [calibration_file,half]=Choose_calibration(calibration_number)

switch calibration_number
    case 1        %%19.9_12
         calibration_file='C:\Users\owner\Documents\Shahar_Neri\friction_videos\19.9\12\C12p.csv';
         half=5645;
    case 2        %%10.10_15a 
         calibration_file='C:\Users\owner\Documents\Shahar_Neri\friction_videos\10.10\15\C15a.csv';
         half=4947;
    case 3        %%10.10_15b 
         calibration_file='C:\Users\owner\Documents\Shahar_Neri\friction_videos\10.10\15\C15b.csv';
         half=4945;
    case 4        %%27.10_17
         calibration_file='C:\Users\owner\Documents\Shahar_Neri\friction_videos\27.10\17\C17.csv';
         half=5434;
    case 5        %%27.10_17p
         calibration_file='C:\Users\owner\Documents\Shahar_Neri\friction_videos\27.10\17\C17p.csv';
         half=6000;
    case 6        %%27.10_20p
         calibration_file='C:\Users\owner\Documents\Shahar_Neri\friction_videos\27.10\20\C20p.csv';
         half=5400;
end


%%% In order to put in a new calibration: 

%%% 1. Open another 'case' in the function
%%% 2. Enter location of calibration file
%%% 3. Insert the middle point 
%%% 4. Update it in the main code 