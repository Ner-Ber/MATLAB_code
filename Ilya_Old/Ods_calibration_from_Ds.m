function [V,d]=Ods_calibration_from_Ds(Ods,ds,baseFig,gtitle)
%[V,d]=Ods_calibration_from_Ds(Ods,ds,baseFig,gtitle)

%The function calibrates Ods according to ds (when the calibartion move is with steps)
% Ods - channel with Oded`s tech. needed calibartion
% ds - displacemnt sensor that calibrats 
%change ds_tresh for better peak detection.
%--------Then
%1.use csaps to interpolate
%2.save the output parameters (V,d,p) at exp_dir\cal\Odscal_ch_094_ex5 (for external ch5) 
%3. use Odsv_2mu to convert involts data to mu

smt=171;%the data is smoothed
minDist=1000; %minimum ditance betwean direvative peaks
ds_tresh=2.3e-4; %treshold for direvative peaks
ds_tresh_max=8.5*ds_tresh; %treshold for max direvative peaks, for the bounderies 

if nargin<4
    gtitle=' ';
end

%-----Cut the better part of the calibration data.  
xx=logical((Ods<4.6).*(Ods>0.4));
ds=ds(xx);
Ods=Ods(xx);

%------smooth the data
ds=my_smooth(ds,smt); %the real displacment (assumed to be converted to microns at 'acq132_convert_raw_to_data'
ds=ds-mean(ds(1:40));
Ods=my_smooth(Ods,smt); %Ods. 

%-----find maxima
abs_diff_ds=abs(diff(Ods)); 
ds_locs=localMaximum(abs_diff_ds,minDist);
ds_locs_tresh=abs_diff_ds(ds_locs)>ds_tresh;
ds_locs_tresh=ds_locs_tresh.*(abs_diff_ds(ds_locs)<ds_tresh_max);
ds_locs_tresh=logical(ds_locs_tresh);
ds_locs=ds_locs(ds_locs_tresh);

%-----make v(d) plot- calibration plot 
sample_index=floor(ds_locs(1:end-1)+1/2*diff(ds_locs));
V=Ods(sample_index); %Voltag
d=ds(sample_index); %displacement

%-----plot 
set(0,'DefaultFigureWindowStyle','docked')

figure(baseFig);
baseFig=baseFig+1;
plot(abs_diff_ds,'.-');
hold all
plot(ds_locs,abs_diff_ds(ds_locs),'o');
hold off
title('Diff(ds)')

figure(baseFig);
baseFig=baseFig+1;
plot(ds,'.-');
hold all
plot(sample_index,ds(sample_index),'o');
hold off
title('ds')

figure(baseFig);
baseFig=baseFig+1;
plot(Ods,'.-');
hold all
plot(sample_index,Ods(sample_index),'o');
hold off
title('Ods')


figure(baseFig);
baseFig=baseFig+1;
% plot(ds,ch5,'.-');
% hold all
%plot(d,V,'o');
plot(V,d,'-o');

title(['calibration plot' gtitle])

