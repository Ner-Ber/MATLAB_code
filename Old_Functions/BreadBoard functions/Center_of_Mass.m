function [X_c,Y_c]=Center_of_Mass(A,InterpFactor) 
B1=sum(A,1);
B2=sum(A,2);

X_c=0;
for i=1:length(B1)
    X_c = X_c+B1(i)*i/InterpFactor;
end
X_c=X_c/sum(B1);

Y_c=0;
for i=1:length(B2)
    Y_c = Y_c+B2(i)*i/InterpFactor;
end
Y_c=Y_c/sum(B2);
    