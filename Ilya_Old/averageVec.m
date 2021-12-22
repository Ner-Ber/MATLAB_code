function [x_out,v_out]=averageVec(x,v)
%Exmple for use: Avergae profiles of contact area . all the profiles are x=0.  
%Matrices should be column orriented


%Find the overlapping x axis
[~,index]=max(min(x));
x_out=x(:,index);
[max_x,~]=min(max(x));
[~,index]=min(abs(x_out-max_x));

if (mean(diff(x_out)) >0)
x_out=x_out(1:index);
else
    x_out=x_out(index:end);
end

%Average the matrice v
[~,index1]=min(abs(x-x_out(1)));
[~,index2]=min(abs(x-x_out(end)));
v_out=zeros(length(x_out),1);
for j=1:length(v(1,:))
v_out=v_out+v(index1(j):index2(j),j);
end
v_out=v_out/length(v(1,:));