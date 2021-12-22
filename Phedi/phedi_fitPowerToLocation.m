function [fitPowerStruct] = phedi_fitPowerToLocation(t,PhediLoc,varargin)
%[fitPowerStruct] = phedi_fitPowerToLocation(t,PhediLoc,t_min, t_max, Shifts)
%
%
% DAEFAULTS
% t_min = 1e-5;
% t_max = 9e-5;
% Shifts = linspace(-2*1e-3,2*1e-3,150);
% selectedPhedis = 'all';


%% set defaults

warning('off');
%% set data
%-- parameters:
t_min_def = 1e-5;
t_max_def = 9e-5;
Shifts_def = linspace(-2*1e-3,2*1e-3,150);
[t_min, t_max, Shifts] = setDefaults4function(varargin, t_min_def,t_max_def, Shifts_def);


%% fit shifted phedis on log scale
%--- create model for mean measurement
MDLs = nan(2,length(Shifts)); %{};
Rsq_vec = nan(1,length(Shifts));
for j = 1:length(Shifts)
    LOGIC_A = t(:)>t_min & t(:)<t_max;
    LOGIC = (t(:)+Shifts(j))>0 & LOGIC_A;
    %--- fit linear
    XX = real(log10(t(LOGIC)+Shifts(j)));
    YY = real(log10(PhediLoc(LOGIC)));
    
    [P,~] = polyfit(XX,YY,1);
    %         slp = (YY(end)-YY(1))/(XX(end)-XX(1));
    %         intr = -slp*XX(end)+YY(end);
    %         P = [slp intr];
    
    MDLs(:,j) = P(:);
    Rsq_vec(j) = r_square(YY, polyval(P,XX));
end

%% save structure
fitPowerStruct = struct('MDLs', MDLs,...
                        'Shifts', Shifts,...
                        'Rsq_vec', Rsq_vec,...
                        't_min', t_min,...
                        't_max', t_max);

warning('on');
end