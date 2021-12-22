function [front]=calc_frontXT_from_XA(phE,tStart,tEnd,smtX,xStart,xEnd)
%The function finds the location in each time step where A drops under the treshold.
%Xstart,Xend in mm important that the data is step like function. Corners
%of the block should be excluded

x=phE.x;
t=phE.t;
A=phE.lines;
dx=mean(diff(x(200:600)));

%--norlmize A and cut the data
firstLine=repmat(mean(A(1:30,:),1),length(A(:,1)),1);
A=A./firstLine;
Treshold=1-5*std(A(1,100:end-100));%2%8

[~,pixStart]=min(abs(phE.x-xStart));
[~,pixEnd]=min(abs(phE.x-xEnd));
[~,t_indexStart]=min(abs(phE.t-tStart));
[~,t_indexEnd]=min(abs(phE.t-tEnd));

x=x(pixStart:pixEnd);
t=t(t_indexStart:t_indexEnd);
A=A(t_indexStart:t_indexEnd,pixStart:pixEnd);

[A,x]=my_smooth(A',smtX,x);
A=A';

%-----create locations for test 
t_test = linspace(tStart,tEnd,10);
for j=1:length(t_test) [~,index_test(j)]=min(abs(t-t_test(j))); end


k=0;
for j=1:length(t)
    
   if(j==1)
            figure;
   end
    %----find x where A crosses threshold
    
    tresh_loc_tmp=find(A(j,:)>0.95 & A(j,:)<Treshold);
    
    if ~isempty(tresh_loc_tmp)
       
        k=k+1;
        %-------Front location by raw data
        tresh_loc(k)=tresh_loc_tmp(end);
        t_index(k)=j;
        
        stage1.t(k)=t(j);
        stage1.A(k)=A(j,tresh_loc(k));
        stage1.x(k)=x(tresh_loc(k));
        
         %---------plot for check
               
        if(~all(index_test-t_index(k))) 
            
            plot(x,A(t_index(k),:),'.-');
            hold all
            plot(stage1.x(k),stage1.A(k),'o');
            %title(num2str(j));
        end
        
        %-----Front location after interpolation of A(:,j)
        %interpolate A(:,t) around the drop and then find exact time when A
        %drops to "Treshold_intrp"
        
%         t_intrp=t(A_loc-20):0.2E-6:t(A_loc+20);
%         A_intrp=csaps(t(A_loc-20:A_loc+20),A(A_loc-20:A_loc+20,j),0.999999999,t_intrp);
%         [~,A_intrp_loc]=min(abs(A_intrp-Treshold_intrp));
%         front_intrp.A(k)=A_intrp(A_intrp_loc);
%         front_intrp.t(k)=t_intrp(A_intrp_loc);
%         front_intrp.x(k)=x(j);
                
    end
    
end



front=stage1;


% front.tRaw=stage2.t;
% front.xRaw=stage2.x;
% front.diff_A_PeakRaw=stage2.diff_A_Peak;
% %---Exclude locations with large derivitive
%
% vRaw=diff(front.xRaw)./diff(front.tRaw);
% include_index=vRaw>-1000;
% include_index=include_index.*(vRaw<3500);
% include_index=logical(include_index);
% include_index=[include_index(1) include_index];
% front.t=front.tRaw(include_index);
% front.x=front.xRaw(include_index);
% front.diff_A_Peak=front.diff_A_PeakRaw(include_index);





