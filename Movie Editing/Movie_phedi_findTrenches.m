function trenchesCoordinates = Movie_phedi_findTrenches(signal,varargin)

% varargin contains the fraction of the whole signal with which ti smooth
% the signal and the shift downwards to use when looking for phedis
%
% updated 11/11/2018:
% trenches are now found by intersection of the above/below zero logical
% and its intersection with the signal-mean, rather than woth the passing
% from negative to positive pixle. 


%% set dafaults:
%--- parameters cell array contains (by this order) the arguments: N, shiftBy
[N_Def, shiftBy_Def] =...
    deal(0.07, 0.5);
Parameters = {N_Def, shiftBy_Def};
for i = 1:length(varargin)
    if ~isempty(varargin{i})
        Parameters{i} = varargin{i};
    end
end
[N, shiftBy] =...
    deal(Parameters{1},Parameters{2});

%%
smoothFactor = round(length(signal)*N)+1;
try
    SmoothedSignal = smooth(signal,smoothFactor);
catch
    disp('');
end
SmoothedSignal = SmoothedSignal-(mean(SmoothedSignal)*shiftBy);

trenchesLogical = SmoothedSignal<0;
%-- old method:
% Dtrenches = diff(trenchesLogical);
% trenchesCoordinates = find([0;Dtrenches(:)]);
%-- updated methos 11/11/2018:
trenchesSign = (trenchesLogical*2-1)*2*max(SmoothedSignal);
XX = 1:length(signal);
[X0,Y0] = intersections(XX,trenchesSign,XX,SmoothedSignal);
trenchesCoordinates = round(X0);

end