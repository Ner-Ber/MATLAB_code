function [front front_intrp]=calc_frontXT_from_TA(phE,smtT,Xstart,Xend)
%The function finds the time at each pixel when normlized A crosses the treshold. To improve the discritization in time, it interpolates with csaps.   
%Xstart,Xend in mm important that the data is step like function. Corners
%of the block should be excluded

Treshold_intrp=0.95;

x=phE.x;
t=phE.t;
lines=phE.lines;

firstLine=repmat(mean(lines(1:30,:),1),length(lines(:,1)),1);
A=lines./firstLine;

[~,pixStart]=min(abs(phE.x-Xstart));
[~,pixEnd]=min(abs(phE.x-Xend));
%pixStart=1;
%pixEnd=length(A(1,:));

k=0;

for j=pixStart:pixEnd
    
    %----find first time A crosses threshold
    A(:,j)=smooth(A(:,j),smtT);
    Treshold=1-8*std(A(1:200,j));
    A_loc=find(A(:,j)>0.9 & A(:,j)<Treshold);
    
    if ~isempty(A_loc)
       
        k=k+1;
        %-------Front location by raw data
        A_loc=A_loc(1);
        front.t(k)=t(A_loc);
        front.A(k)=A(A_loc,j);
        front.x(k)=x(j);
        
        %-----Front location after interpolation of A(:,j)
        %interpolate A(:,t) around the drop and then find exact time when A
        %drops to "Treshold_intrp"
        
        t_intrp=t(A_loc-20):0.2E-6:t(A_loc+20);
        A_intrp=csaps(t(A_loc-20:A_loc+20),A(A_loc-20:A_loc+20,j),0.999999999,t_intrp);
        [~,A_intrp_loc]=min(abs(A_intrp-Treshold_intrp));
        front_intrp.A(k)=A_intrp(A_intrp_loc);
        front_intrp.t(k)=t_intrp(A_intrp_loc);
        front_intrp.x(k)=x(j);
                
    end
end


x_test=40:10:180;
for j=1:10 [~,index(j)]=min(abs(x-x_test(j))); end

figure;
plot(t,A(:,index),'.-');
hold all;

for j=1:10 [~,index(j)]=min(abs(front_intrp.x-x_test(j))); end
plot(front.t(index),front.A(index),'o');
plot(front_intrp.t(index),front_intrp.A(index),'o');


