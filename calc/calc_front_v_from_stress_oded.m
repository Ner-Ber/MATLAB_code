function [x trise tmax tdown]=acq132_event_front_v_from_stress(exp_dir,event,from,to,smt,baseFig)
%poisson ratio
sigma=1/3;
%youngs modulus
E=3E+9;

if (nargin<6)
    baseFig=1;
end

%------------
sg_093_order=[15:-1:1 31:-1:17];
sg_094_order=15:-1:1;

NFactor= 2267.96185 / (0.003 * 15.031 * 84.3333);
FFactor= 1000 / (0.003 * 15.031 * 84.3333);
%!!!!!!!!!!!!!!!! the sample is supposedly 15cm long - need to check this more accurately !!!!!!!!!!!!!!!!
x_sg=150-[6.0421   14.8074   23.1472   31.1466   38.6354   47.4858   56.0809   71.5691   85.2702   94.2057  103.5667  119.2251  126.6288  136.4153  144.67];
%--------------

path_093=['C:\Frics\' exp_dir '\acq132_093\multivent\'];
path_094=['C:\Frics\' exp_dir '\acq132_094\multivent\'];
%---------------read event data and offset correction
host='acq132_093';
[t,ch_093]=acq132_event_read(host,path_093,1);%using event1 for offset correction
ch_offset=mean(ch_093);
[t,ch_093]=acq132_event_read(host,path_093,event);
ch_093=ch_093-repmat(ch_offset,length(ch_093(:,1)),1);%offset correction// comment this line for disable

host='acq132_094';
[t,ch_094]=acq132_event_read(host,path_094,1);%using event1 for offset correction
ch_094(:,end-2:end)=-1*ch_094(:,end-2:end);
ch_offset=mean(ch_094);
[t,ch_094]=acq132_event_read(host,path_094,event);
ch_094(:,end-2:end)=-1*ch_094(:,end-2:end);
ch_094=ch_094-repmat(ch_offset,length(ch_094(:,1)),1);%offset correction// comment this line for disable

%--------------------orgenize the data and convert units
sg=[ch_093(:,sg_093_order) ch_094(:,sg_094_order)]; %ordered strain gages

clear ch_094 ch_093;

ExcV=0.6;
Gain=470;
%--------------------calculate strain stress
Uxy=-(sg(:,1:3:end)-sg(:,3:3:end))/2/Gain/ExcV;%calculated shear strain
Uyy=sg(:,2:3:end)/Gain/ExcV;%calculated normal strain
Uxx=(sg(:,1:3:end)+sg(:,3:3:end)-Uyy)/Gain/ExcV;%calculated transvers strain

%now calc stresses using hook law and assuming Szz=0
Uzz=-sigma/(1-sigma)*(Uxx+Uyy);
Sxx=E/((1+sigma)*(1-2*sigma))*((1-sigma)*Uxx+sigma*(Uyy+Uzz));
Syy=E/((1+sigma)*(1-2*sigma))*((1-sigma)*Uyy+sigma*(Uxx+Uzz));
Sxy=E/(1+sigma)*Uxy;

sSxy=smoothts((-Sxy./Syy)','b',smt)';
firstLine=repmat(mean(sSxy(1:from,:)),length(t),1);

% sSxy=abs(sSxy-firstLine);
sSxy=(sSxy-firstLine);

%thresholds
rmSxy=2*std(Sxy(1+smt:from-smt-1,:));
rmsLine=repmat(rmSxy,length(t),1);

% meanSxy=mean(sSxy(1+smt:from-smt-1,:));
% meanLine=repmat(meanSxy,length(t),1);

% threshLine=meanLine+rmsLine;

%--------------------plot
set(0,'DefaultFigureWindowStyle','docked')
% hold all
% figure(baseFig);

for i=1:length(sSxy(1,:))
%     plot(from:to,sSxy(from:to,i),'.-');
%     hold all
    %find cuts
    iup=find(sSxy(from:to,i)>rmsLine(from:to,i),1,'first')+from;
    idown=find(sSxy(from:to,i)<-rmsLine(from:to,i),1,'first')+from;
        
%     plot(from:to,rmsLine(from:to,i),'-');
%     plot(from:to,-rmsLine(from:to,i),'-');
    if (~isempty(iup))
%         plot([iup iup],[min(sSxy(from:to,i)) max(sSxy(from:to,i))],'-');
        trise(i)=iup;
    else
        trise(i)=NaN;
    end
    if (~isempty(idown))
%         plot([idown idown],[min(sSxy(from:to,i)) max(sSxy(from:to,i))],'-');
        tdown(i)=idown;
    else
        tdown(i)=NaN;
    end
    
    %find middle peak
    if (~isempty(idown))
        if (~isempty(iup) && tdown(i)>trise(i))
            tmax(i)=find(sSxy(trise(i):tdown(i),i)==max(sSxy(trise(i):tdown(i),i)),1,'first')+trise(i)-1;
        else
            tmax(i)=find(sSxy(from:tdown(i),i)==max(sSxy(from:tdown(i),i)),1,'first')+from-1;
        end
        
%          plot([tmax(i) tmax(i)],[min(sSxy(from:to,i)) max(sSxy(from:to,i))],'-');
    else
        tmax(i)=NaN;
    end
       
%     if (~isempty(iup))
%         if (~isempty(idown))
%         end
%     end
    
%     title(sprintf('Sg x - %f',x_sg(i)));
%     hold off
%     if (i<length(sSxy(1,:)))
%         waitforbuttonpress
%     end
    
end
% title('Shear Stress');

x=x_sg;

%get rid of obviously wrong points
trise(x==min(x))=NaN;

figure(baseFig);
hold off

plot(trise,x,'.','DisplayName','trise');
hold all
plot(tmax,x,'.','DisplayName','tmax');
plot(tdown,x,'.','DisplayName','down');
hold off
xlabel('t (\musec)');
xlim([from to]);
ylabel('x (mm)');
title(sprintf('Event# - %d, smt - %d',event,smt));
legend;

figure(baseFig+1);
hold off


% tbmax=interp1(x(1:end-3),tmax(1:end-3),x-0.1,'spline');
tamax=interp1(x(1:end-3),tmax(1:end-3),x+1,'spline');
vp=1./(tamax-tmax);

for i=1:length(vp)
  tx=vp(i)*((1:length(t))-tmax(i));
  plot(-tx(from:to),sSxy(from:to,i)/sSxy(tmax(i),i),'.-');
  hold all
end
xlim([-50 250])
ylim([-0.6 1.1])
hold off

figure(baseFig+2);
hold off

% tbdown=interp1(x(1:end-3),tdown(1:end-3),x-1,'spline');
tadown=interp1(x(1:end-3),tdown(1:end-3),x+1,'spline');
vd=1./(tadown-tdown);

for i=1:length(vd)
  tx=vd(i)*((1:length(t))-tdown(i));
  plot(-tx(from:to),sSxy(from:to,i),'.-');
  hold all%/sSxy(tdown(i),i)
end
xlim([-50 250])
% ylim([-0.6 1.1])
hold off

figure(baseFig+3);
plot(x,vp*1000,'.');

figure(baseFig+4);
hold off
for i=1:length(x)
plot((1:length(t))-tmax(i),sSxy(:,i),'.-');
hold all
end
hold off



