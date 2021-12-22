function plot_Sxy_Ds_A_normlized(exp_dir,event,sg_num,ch_094_ex_num,baseFig)
%The function calculates the Stresses and displacement at perticular sg.
%And normolize the data. 
Sxy_smt=1;
sensF=-2.1;%[kg/micron]

acq132_event_get_data(exp_dir,event,'max','max',1,'Sxy','ch_094_ex','F');

Ds=ch_094_ex(:,ch_094_ex_num); %displacment in lab frame 
slip=Ds-(F-mean(F(1:20)))/sensF; %slip relative to the stage
Sxy=smooth(Sxy(:,sg_num),Sxy_smt);

%--find first drop index
smt=51;
Ds_abs=abs(smooth(Ds,smt));
min_i=localMaximum(Ds_abs,1600);
min_i_i=find(min_i>smt);   %smooth is not good at the bounderies
min_i=min_i(min_i_i);
min_i_tresh=find(Ds_abs(min_i)>1); %treshold for the maximum
min_i=min_i(min_i_tresh);
min_i=min_i(1);

%--normolize stress and slip by the value of the first drop
Ds_n=Ds/abs(mean(Ds(min_i-10:min_i+10))); 
slip_n=slip/abs(mean(slip(min_i-10:min_i+10)));
Sxy_n=norm_data(Sxy);
Sxy_n=Sxy_n/abs(mean(Sxy_n(min_i-10:min_i+10)));

%----plot the data

figure(baseFig);
baseFig=baseFig+1;
plot(slip,Sxy,'.-');
title('Sxy(slip)');

figure(baseFig);
baseFig=baseFig+1;
plotyy(t,Sxy,t,slip,'plot');
%to get axes from plotyy use - ax=findobj(gcf,'Type','axes');


figure(baseFig);
baseFig=baseFig+1;
plot(t,Sxy_n,'.-');

% [x_sg_spl,Fx,Syy_spl,Sxy_spl,Sxx_spl,dSxx_spl]=calc_Fx_stress_splined(x_sg,Sxx,Syy,Sxy);
% sg_num=find(x_sg_spl==x_sg(sg_num));
% Fx=norm_data(Fx(:,sg_num));
% Fx=Fx/abs(mean(Fx(I-10:I+10)));
% plot(t,Fx,'.-');

hold all
plot(t,[slip_n Ds_n],'.-');
title(['ch094ex num - ' num2str(ch_094_ex_num)]);
hold off


