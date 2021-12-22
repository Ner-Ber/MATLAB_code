function [vSmt]=my_smooth2(x,v,smtX)

%The function smoothes vector v over range of smtX according to x vector.
%good for ruptures that propagate very slowly and then raidly accelerat. 

for j=1:length(x)
    index1=find(x-x(j)<-smtX/2,1,'last');
    if(isempty(index1))
        index1=1;
    end
    index2=find(x-x(j)>smtX/2,1,'first');
    
    if(isempty(index2))
        index2=length(x);
    end
    vSmt(j)=mean(v(index1:index2));
end


