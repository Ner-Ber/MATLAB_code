%% lad figure and data and get axes
[Cd, Cs, Cr, nu, ro, E, mu, ~, PlaneStrain, tau_p, Xc0] = CrackSolutionMaterialProperties;
openfig('C:\Users\NB\Google Drive\JAY_lab\WRITING\FIGURES\fig4_tau\subplots\bothHighAndLowCrackTip\Vx_and_A_high.fig');
F = gcf;
A = F.Children.Children(arrayfun(@(A) strcmpi('mean phedi',A.DisplayName),F.Children.Children));
x_meanPhedi = A.XData;
y_meanPhedi = A.YData;
delete(F);      % close open figure
%% replot
explain_tau_nl_F = figure;
ax = axes;
VxLine = plot(x_meanPhedi, y_meanPhedi);

%% load data of points
A1 = load('F:\NeriFricsData\Frics\2018-10-10\DAandREsidualsStruct.mat');
myStruct20181010 = A1.myStruct;
%--- add phedi data:
% phediStruct = load('F:\NeriFricsData\Frics\2018-10-10\allPhediStructures_17186.mat');
% phediStruct = phediStruct.PhediStructCell;
%--- get specific points and data:
tau_nl = myStruct20181010(7).yieldSxy2ndAll(myStruct20181010(7).Events==12);
X_nl = myStruct20181010(7).rise2XAll(myStruct20181010(7).Events==12);
Gamma = myStruct20181010(7).Gamma(myStruct20181010(7).Events==12);
Cf = myStruct20181010(7).Cf_vec(myStruct20181010(7).Events==12);

%% --- create the LEFM model
[cohesiveFun,PlaneFlag] = deal('exp','PlaneStrain');
% Nx = 30;
% X_vec = linspace(-2,0,Nx);
% DX = X_vec(2)-X_vec(1);
v = Cf/Cr;
L = 10*3.2e-3;
interfaceEpsilon = 1e-9;
X =  -(-0.1:1e-7:0.1);
solAtInter = CrackSolutionForh_GammaChange(v, -0.5, interfaceEpsilon, Gamma, X);
solAtInter.Sxy(solAtInter.x<=0) = 0;
yyaxis right
SxyLine = plot(solAtInter.x,solAtInter.Sxy);

%-- -set proportions:
ylim([0,13e5]);
xlim([-1,1]*1e-1)
ax.YAxis(2).Limits(1) = ax.YAxis(2).Limits(2)*ax.YAxis(1).Limits(1)/ax.YAxis(1).Limits(2);
ax.XAxisLocation='origin';

%--- ad measurement markings:
% yyaxis left;
% hold on;
% plot(X_nl,0,'o');
yyaxis right;
hold on;
plot([1,1]*X_nl,[0,tau_nl],'-');






