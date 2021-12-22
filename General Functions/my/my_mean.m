function [tt,yy]= my_mean(t,y)
%example for use: you have normolized self similar signal over the different sg, and
%already you applied the offset.Now you want to get the mean signal. 
%The function cuts the data appropriatly and returns the mean signal

tStart=max(t(1,:));
tEnd=min(t(end,:));

yy=0;
tt=0;
for j=1:length(y(1,:))
    [~,index1]=min(abs(t(:,j)-tStart));
    [~,index2]=min(abs(t(:,j)-tEnd));
    yy=yy+y(index1:index2,j);
    tt=tt+t(index1:index2,j);
end
yy=yy/length(y(1,:));
tt=tt/length(y(1,:));
