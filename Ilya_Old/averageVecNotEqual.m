function [x_out,v_out]=averageVecNotEqual(x,v)
%Average vectors with different coordinates.Example: different strain gages
%after going from t to x.
%Matrix x,v should be column oriented.


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

for j=1:length(v(1,:))
vi(:,j)=interp1(x(:,j),v(:,j),x_out,'spline');
end
%Average the matrice v
v_out=mean(vi,2);