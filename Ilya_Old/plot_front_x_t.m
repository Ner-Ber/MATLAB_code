function plot_front_x_t(x_sg,baseFig) 

f=calc_front_x_t(x_sg); %front struct
v=diff(f.x)./diff(f.t);%[m/sec]
v_x=f.x(1:end-1)+diff(f.x)/2;% the mean coordiante of velocity calculation
v_error=1./diff(f.t).*(0.085^2+(0.001*v_x).^2).^(1/2)/2; %1 pix error->0.085 mm in space, 0.001 msec error in time

%---interp the mean velocity
v_spl_x=f.x(2:end-1);
v_spl=interp1(v_x,v,v_spl_x,'linear'); 
v_spl_x=[f.x(1) v_spl_x f.x(end)]; %Don't want to waste last point. So , approx velocity at the last point (no interpolation) 
v_spl=[v(1) v_spl v(end)];

%---norm the data
[t,Sxy]=my_get_axis;
Sxy_error=std(Sxy(1:100,f.sg_num))/2;
Sxy_n=(Sxy(:,f.sg_num)-repmat(f.Sxy_i,length(t),1))./repmat(f.Sxy_p-f.Sxy_i,length(t),1);
x_n=(repmat(t,1,length(f.t))-repmat(f.t,length(t),1)).*repmat(v_spl,length(t),1);%[m/s]*[msec]->[mm] 

%--plot data
figure(baseFig);
baseFig=baseFig+1;
plot(x_n,Sxy_n,'.-')

figure(baseFig);
baseFig=baseFig+1;
errorbar(v_x,v,v_error,'.-')
hold all
plot(v_spl_x,v_spl,'o')
hold off

figure(baseFig);
baseFig=baseFig+1;
errorbar(f.x,f.Sxy_p-f.Sxy_i,Sxy_error,'.-')
