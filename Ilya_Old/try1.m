xc=0.5E-3;%[m]
v=4000;
dt=10^-6;

a=xc/(v*dt);
x=(-7:1/20:70);

A=tanh(x);

A_dt=a*log(cosh(x)./cosh(x-a^-1));

%figure;

plot(x*xc,A,'.-');
hold all;
%plot((x+a^-1)*xc,A,'.-')
plot(x*xc,A_dt,'.-');

% figure;
% 
% plot(x(2:end),diff(A),'.-');
% hold all;
% plot(x(2:end)+a^-1,diff(A),'.-')
% plot(x(2:end),diff(A_dt),'.-');