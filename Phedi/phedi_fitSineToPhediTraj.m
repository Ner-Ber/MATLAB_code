function myfitStruct = phedi_fitSineToPhediTraj(PhediData,varargin)
%myfitStruct = phedi_fitSineToPhediTraj(PhediData,preKinkTime,numOfSin,DoPlot)
%
% only mandatory input is 'PhediData'.
% phedi_fitSineToPhediTraj will take the location of [hedis pre the
% 'preKinkTime' and fit two sins (default) or a single sin to it.
%
% DEFAULTS:
% preKinkTime = -3e-4
% numOfSin = 2
% DoPlot = 0

%% set defaults
[preKinkTime,numOfSin,DoPlot] = setDefaults4function(varargin,-0.025/PhediData.Cf,2,0);
%% set the data
t_mins_t_tips = bsxfun(@minus,PhediData.timeVec(:),PhediData.t_tips(:)');  % t-t_tip
PhediLoc = PhediData.PhediLocation;     % phedilocation vectors (columns of matrix)
preLogical = t_mins_t_tips<preKinkTime; % locations where is relevant to fit sine (logical)
[L, N] = size(t_mins_t_tips);              % number of phedis
pre_t_mins_t_tips_Cell = mat2cell(t_mins_t_tips,L,ones(1,N));
prePhediLocCell = mat2cell(PhediLoc,L,ones(1,N));
preLogicCell = mat2cell(preLogical,L,ones(1,N));

preTimeCellWithNan = cellfun(@(A,B) A(B), pre_t_mins_t_tips_Cell, preLogicCell, 'UniformOutput',0);    % time vectors for fitting
preLocCellWithNan = cellfun(@(A,B) A(B), prePhediLocCell, preLogicCell, 'UniformOutput',0);            % location vectors for fitting
preLocCell = cellfun(@(A) A(~isnan(A)), preLocCellWithNan, 'UniformOutput',0);
preTimeCell = cellfun(@(A,B) A(~isnan(B)), preTimeCellWithNan, preLocCellWithNan, 'UniformOutput',0);

%% find amplitude and shift
Amps = cellfun(@(A) peak2peak(A,1),preLocCell);
Shifts = cellfun(@(A) mean(A,1,'omitnan'),preLocCell);

%% fit per vector
options = optimset('MaxFunEvals',1e5,'Display','off');

myfitCell = cell(1,N);
Phases = nan(1,N);
%--- prepare the figure if needed
if DoPlot
    n = ceil(3*N/12);
    m = ceil(N/n);
    figure;
    ha = tight_subplot(n,m,[.02 .03],[.1 .01],[.01 .01]);
end
for i = 1:N
    t = preTimeCell{i};
    y = preLocCell{i};
    [~,min_y_Idx] = min(y,[],'omitnan');
    [~,max_y_Idx] = max(y,[],'omitnan');
    phase_i = t(min_y_Idx);
    Shift_i = Shifts(i);
    Amp_i = Amps(i);
    initial_omega = pi/abs(t(min_y_Idx)-t(max_y_Idx));
    switch numOfSin
        case 1
            modelfun = @(b,x) Shift_i-0.5*abs(b(3))*cos(b(1)*(x - b(2)));
            fcn = @(B) sum((modelfun(B,t) - y).^2);
            b0 = [pi/abs(t(min_y_Idx)-t(max_y_Idx)), phase_i, Amp_i];
            
            myfitCell{i} = fminsearch(fcn, b0, options);
            
            if DoPlot
                axes(ha(i));
                plot(t,y,'.');
                hold on;
                plot(t,modelfun(myfitCell{i},t),'r','LineWidth',1.5);
                title(['i=',num2str(i)],'FontSize',3);
            end
            
        case 2
            Sin1_model = @(b,x) Shift_i-0.5*abs(b(3))*cos(b(1)*(x - b(2)));
            fcn_sin1 = @(B) sum((Sin1_model(B,t) - y).^2);
            b01 = [pi/abs(t(min_y_Idx)-t(max_y_Idx)), phase_i, Amp_i];
            
            s1 = fminsearch(fcn_sin1, b01, options);
            y2 = y - Sin1_model(s1,t);
            
            Sin2_model = @(b,x) b(1)+b(2)*cos(b(3)*(x - b(4)));
            fcn_sin2 = @(B) sum((Sin2_model(B,t) - y2).^2);
            b02 = [0, Amp_i/5, initial_omega*4, 0];
            s2 = fminsearch(fcn_sin2, b02, options);
            
            totalModel = @(b,x) Shift_i-0.5*abs(b(3))*cos(b(1)*(x - b(2)))  +...
                b(4)+b(5)*cos(b(6)*(x - b(7)));
            myfitCell{i} = [s1(:)', s2(:)'];
            
            if DoPlot
                axes(ha(i));
                plot(t,y,'.');
                hold on;
                plot(t,totalModel(myfitCell{i},t),'r','LineWidth',1.5);
                title(['i=',num2str(i)],'FontSize',3);
            end
            
        otherwise
            error('number of sines invalid');
    end
    
    %     myfitCell{i} = fminsearch(fcn, b0, options);
    Phases(i) = phase_i;
    
    
end

%% save stuff
myfitStruct = struct;
myfitStruct.b_param = myfitCell;
myfitStruct.Phases = Phases;
myfitStruct.Shifts = Shifts;
myfitStruct.Amps = Amps;
switch numOfSin
    case 1
        myfitStruct.model= @(b,shift,x) shift-0.5*abs(b(3))*cos(b(1)*(x - b(2)));
    case 2
        myfitStruct.model= @(b,shift,x) shift-0.5*abs(b(3))*cos(b(1)*(x - b(2)))  +...
            b(4)+b(5)*cos(b(6)*(x - b(7)));
end



end