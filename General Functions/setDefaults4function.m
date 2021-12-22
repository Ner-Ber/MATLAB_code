function [varargout] = setDefaults4function(vararginInput,varargin)
% [variableNames] = setDefaults4function(varargin,VariableDefaualtValues)
%
%setDefaults4function will set dafault for a mother function.
% vararginInput - is the varargin of the mother function.
% varargin - will contain the default values wanted (all of them!)
% varargout - will be the arguments in the mother function which you are
% setting the defaults for.
% NUMBER OF OUTPUTS MUST EQUAL length(varargin) !!!

varargout = cell(size(varargin));

if length(vararginInput)==1
    if iscell(vararginInput{1})
        vararginInput = vararginInput{1};
    end
end

if ~isempty(vararginInput)
    for i = 1:length(vararginInput)
        if ~isempty(vararginInput{i})
            varargout{i} = vararginInput{i};
        else
            varargout{i} = varargin{i};
        end
    end
end
varargout(length(vararginInput)+1:end) = varargin(length(vararginInput)+1:end);
end