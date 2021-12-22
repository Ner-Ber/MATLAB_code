function [leftLim rightLim]=vdsCalcLims(dirname,interval)
eval(sprintf('cd %s',dirname))
ims=vdsReadImsSlow(interval);
 cd ..

 lines1=[ims.Line1]';
left=zeros(length(lines1(:,1)),1);
right=zeros(length(lines1(:,1)),1);
for i=1:length(lines1(:,1))
    left(i)=find(lines1(i,:)>0,1,'first');
    right(i)=find(lines1(i,:)>0,1,'last');
end

leftLim=min(smooth(left,100));
rightLim=max(smooth(right,100));

