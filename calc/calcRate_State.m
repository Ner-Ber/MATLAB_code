function calcAfromSlipV(path)

global pathMfile;

pathMfile='C:\Frics\2013-03-17\19-17-43\e7\try2.mat';

%t=0:10^-7:2*10^-4;
% phi_0=5;%[sec]
% phi_x=10^-4;%[sec]
% d=4*10^-6;%[m]
% V_x=3*d/phi_x;%[m/s]
% beta=0.04;


phi_0=6;%[sec]
phi_x=10^-5;%[sec]
beta=0.05;

load(pathMfile);

t_start=slipV_t(1);
t_end=slipV_t(end);

%--check slipV interpolation
v=calc_v(slipV_t);
figure(3);
plot(slipV_t,v,'.-');

%----solve Ode
[T,phi] = ode45(@phi_evol,[t_start t_end],phi_0);
A=1+beta*log(1+phi/phi_x);
figure(4);
plot(T,A/A(1),'.-');

%---step v solution
% d=0.9E-6;
% V_x=d/phi_x;
% V=0.06;
% t=slipV_t;
% A=1+beta.*log(1+V_x/V+(phi_0/phi_x-V_x/V)*exp(-V*(t-t(1))/d));
% figure(2);
% hold all;
% plot(t,A/A(1),'.-');


function dy = phi_evol(t,phi)
d=0.5*10^-6;%[m]
v=calc_v(t);
dy = 1-phi*v/d;%dy = 1-phi*v/d;

function v=calc_v(t)
global pathMfile
load(pathMfile);
v=csaps(slipV_t,slipV,1,t);
%v=0.06*ones(1,length(t));





