function [front]=calc_frontXT_from_diffT_A(phE,smtT,Xstart,Xend)
%The function finds the time at each pixel when time derivitive of A is maximum.
%Xstart,Xend in mm important that the data is step like function. Corners
%of the block should be excluded


x=phE.x;
t=phE.t;
dt=mean(diff(t));
lines=phE.lines;

%--norlmize A and cut the data
firstLine=repmat(mean(lines(1:30,:),1),length(lines(:,1)),1);
A=lines./firstLine;
[A,t]=my_smooth(A,smtT,t);

[~,pixStart]=min(abs(phE.x-Xstart));
[~,pixEnd]=min(abs(phE.x-Xend));

x=x(pixStart:pixEnd);
A=A(:,pixStart:pixEnd);

%-----create locations for test at stage 2 loop
x_test = (20:10:170);
for j=1:length(x_test) [~,index_test(j)]=min(abs(x-x_test(j))); end


%-----
diff_A=-diff(A);
diff_t=t(2:end)+(t(2)-t(1))/2;


%-----First stage: find  first maxima in the time derivitive

for j=1:length(x)
    treshold=5*std(diff_A(1:100,j));
    max_loc=localMaximum(diff_A(:,j),500); %100 was chosen ampiracly.for example 1 is not working
    above_tresh_loc=diff_A(max_loc,j)>treshold; %returns binary vector
    above_tresh_loc=max_loc(above_tresh_loc);
    diffPeakLoc(j)=above_tresh_loc(1);
    
    stage1.diff_A_Peak(j)=diff_A(diffPeakLoc(j),j);
end
stage1.t=diff_t(diffPeakLoc);
stage1.x=x;



%-------Second Stage:  find the time by parabolic fit
v=diff(stage1.x)./smooth(diff(stage1.t),51)';
v=[v(1) v];

for j=1:length(x)
    
    %deffine the number of point to fit. should depend on front velocity
    %may be imroved
    if (v(j)<250)
        delta=11;
    else
        delta=9;
    end
    
    if (diffPeakLoc(j)-delta <1)||(diffPeakLoc(j)+delta >length(diff_t)) %check if tring to access over the limits
        
        stage2.t(j)=diff_t(diffPeakLoc(j));
        stage2.diff_A_Peak(j)=diff_A(diffPeakLoc(j));
        
    else
        
        for_fit.t=diff_t(diffPeakLoc(j)-delta:diffPeakLoc(j)+delta);
        for_fit.diff_A=diff_A(diffPeakLoc(j)-delta:diffPeakLoc(j)+delta,j);
        p= polyfit(for_fit.t,for_fit.diff_A,2);
        
        stage2.t(j)=-p(2)/2/p(1);
        stage2.diff_A_Peak(j)=polyval(p,stage2.t(j));
        
        %if the results too different. take the raw data
        
        if ( abs(stage2.t(j)-stage1.t(j))/dt) >10
            stage2.t(j)=stage1.t(j);
            stage2.diff_A_Peak(j)=stage1.t(j);
        end
        
        stage2.x(j)=x(j);
        %---------plot for check
        if(j==1) 
            figure; 
            hold all; 
        end
        
        if(~all(index_test-j))
            plot(diff_t,diff_A(:,j),'.-');
           % plot(stage1.t(j),stage1.diff_A_Peak(j),'o');
            plot(for_fit.t,polyval(p,for_fit.t),'.-')
            %plot(stage2.t(j),stage2.diff_A_Peak(j),'o');
            
        end
        %-----------
    end
    
end


%front=stage2;

front.tRaw=stage2.t;
front.xRaw=stage2.x;
front.diff_A_PeakRaw=stage2.diff_A_Peak;
%---my_robust

vRaw=diff(front.xRaw)./diff(front.tRaw);
include_index=vRaw>-1000;
include_index=include_index.*(vRaw<3500);
include_index=logical(include_index);
include_index=[include_index(1) include_index];
front.t=front.tRaw(include_index);
front.x=front.xRaw(include_index);
front.diff_A_Peak=front.diff_A_PeakRaw(include_index);


%--plot the the data and the locations
% x_test = (20:10:170);
% for j=1:length(x_test) [~,index(j)]=min(abs(x-x_test(j))); end
% 
% figure;
% plot(diff_t,diff_A(:,index),'.-');
% hold all;
% plot(front.t(index),front.diff_A_Peak(index),'o');
% 
 legend(num2str(x_test'));legend off;


