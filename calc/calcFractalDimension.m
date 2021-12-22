function [C,r]=calcFractalDimension(v)

%creat logistic map
%lambda=3.5699456;
lambda=3.57;
f=@(x)lambda*x.*(1-x);
x_vec=zeros(1,1e4);
x_vec(1)=0.5;

for j=2:length(x_vec);
    x_vec(j)= f(x_vec(j-1));
end
% ----------------
v=x_vec;
%for 1D
r=logspace(-5,-0,500);
C=r*0;

for k=2:length(v)-1
    vtmp=abs( [v(1:k-1) v(k+1:end)]-v(k) );
    
    for l=1:length(r)
        C(l)=C(l)+sum(vtmp < r(l));
     end
    
end
C=C/length(v)^2;
