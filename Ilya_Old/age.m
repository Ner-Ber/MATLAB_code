function dy = age(t,y)
w=0.1;
b=0.01;
a=0.2;
S=1;
dS=1*S;

dy =1/w-a/b*y*(dS*cos(t)/(S+dS*sin(t)));