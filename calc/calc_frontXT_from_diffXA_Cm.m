function [front]=calc_frontXT_from_diffXA_Cm(phE,tStart,tEnd,smtX,xStart,xEnd)
%The function finds the time at each pixel when time derivitive of A is maximum.
%Xstart,Xend in mm important that the data is step like function. Corners
%of the block should be excluded
%stage1-find the time of the peak at the time derivitive in each pixel.
%about that time calculate the "centar of mass". the domain for the
%calcluation is defined by "delta".
%Sometimes need to change the delta and combine diferent

x=phE.x;
t=phE.t;
A=phE.lines;
dx=mean(diff(x(200:600)));

%--norlmize A and cut the data
firstLine=repmat(mean(A(1:30,:),1),length(A(:,1)),1);
A=A./firstLine;

[~,pixStart]=min(abs(phE.x-xStart));
[~,pixEnd]=min(abs(phE.x-xEnd));
[~,t_indexStart]=min(abs(phE.t-tStart));
[~,t_indexEnd]=min(abs(phE.t-tEnd));

x=x(pixStart:pixEnd);
t=t(t_indexStart:t_indexEnd);
A=A(t_indexStart:t_indexEnd,pixStart:pixEnd);

[A,x]=my_smooth(A',smtX,x);
A=A';


%-----
diff_A=diff(A,1,2);
diff_x=x(2:end)+(x(2)-x(1))/2;


%-----First stage: find  first maxima in the time derivitive
treshold=5*std(diff_A(1,100:end-100));
k=0;

for j=1:length(t)
    
    max_loc=localMaximum(diff_A(j,:),500); %500 was chosen ampiracly.for example 1 is not working
    above_tresh_loc=diff_A(j,max_loc)>treshold; %returns binary vector
    
    if (nnz(above_tresh_loc))
      k=k+1;  
        above_tresh_loc=max_loc(above_tresh_loc);
        diffPeakLoc(k)=above_tresh_loc(end);
        t_index(k)=j;
        
        stage1.diff_A_Peak(k)=diff_A(j,diffPeakLoc(k));
        stage1.t(k)=t(j);
        stage1.x(k)=diff_x(diffPeakLoc(k));
        
      end
    
end

%-----create locations for test 
t_test = linspace(tStart,tEnd,10);
for j=1:length(t_test) [~,index_test(j)]=min(abs(stage1.t-t_test(j))); end


%-------Second Stage:  find the center of mass
% v=diff(stage1.x)./smooth(diff(stage1.t),51)';
% v=[v(1) v];

for j=1:length(stage1.t)
    
    %deffine the number of point to fit. should depend on front velocity
    %may be imroved
    % if (v(j)<250)
    delta=17;%11
    %else
    %   delta=7;
    %end
    
    if (diffPeakLoc(j)-delta <1)||(diffPeakLoc(j)+delta >length(diff_x)) %check if trying to access over the limits
        
        stage2.x(j)=stage1.x(j);
        stage2.t(j)=t(j);
        stage2.diff_A_Peak(j)=stage1.diff_A_Peak(j);
        
        
    else
        
        
        cm_x=diff_x(diffPeakLoc(j)-delta:diffPeakLoc(j)+delta); %cm-Center of mass
        cm_diff_A=diff_A(t_index(j),diffPeakLoc(j)-delta:diffPeakLoc(j)+delta);
        stage2.x(j)=sum(cm_x.*cm_diff_A)/(sum(cm_diff_A));
        stage2.diff_A_Peak(j)=sum(cm_diff_A)/length(cm_diff_A);
        
        %if the results too different. take the raw data
        
        if ( abs(stage2.x(j)-stage1.x(j))/dx) >10
            stage2.x(j)=stage1.x(j);
            stage2.diff_A_Peak(j)=stage1.diff_A_Peak(j);
        end
        
        stage2.t(j)=t(j);
        %---------plot for check
        if(j==1)
            figure;
        end
        
        if(~all(index_test-j))
            
            plot(diff_x,diff_A(t_index(j),:),'.-');
            hold all
            plot(stage1.x(j),stage1.diff_A_Peak(j),'o');
            plot(cm_x,cm_diff_A);
            plot(stage2.x(j),stage2.diff_A_Peak(j),'o');
            title(num2str(j));
        end
        %-----------
    end
    
end


front=stage2;

front.tRaw=stage2.t;
front.xRaw=stage2.x;
front.diff_A_PeakRaw=stage2.diff_A_Peak;
%---Exclude locations with large derivitive

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
%legend(num2str(x_test'));legend off;


