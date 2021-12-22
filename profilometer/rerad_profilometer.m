
% test=dir;
% filename=test(9).name;

test=dir_names;
filename=test{1};


A=importdata(filename);

Altitude=A.data(:,2);  
xx=A.data(:,3)*1e-3;   %[mm]


fixed_xx = xx(Altitude>100);
fixed_Altitued = Altitude(Altitude>100);
fixed_Altitued = detrend(fixed_Altitued);
fixed_Altitued=-fixed_Altitued;

figure; hold all
plot(xx,Altitude,'.')
plot(fixed_xx-fixed_xx(1),fixed_Altitued)

title(strrep(filename(1:20),'_',' '));
xlabel('x [mm]')
ylabel('Altitude [\mum]')