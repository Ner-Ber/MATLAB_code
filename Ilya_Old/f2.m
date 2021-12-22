function ff=f2(C)

%the analog of gamma=(1-(v/cs)^2)^0.5 in antiplane

Cs=1370;
Cd=2700;
nu=1/3;

alpha_s=(1-(C/Cs).^2).^0.5;
alpha_d=(1-(C/Cd).^2).^0.5;
D=4*alpha_s.*alpha_d-(1+alpha_s.^2).^2;

ff=D*(1-nu)*Cs^2./alpha_s./C.^2; %inverse of the f function

