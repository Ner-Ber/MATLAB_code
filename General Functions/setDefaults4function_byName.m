function Param = setDefaults4function_byName(defaultInputs,vararginInput)
% function setDefaults4function_byName(defaultInputNames,varargin)
%
% defaultInputs should contain default inputs in the form of a cell array as such:
% {'inputName_a',defaultInputValue_a;...
% 'inputName_b',defaultInputValue_b;...
% 'inputName_c',defaultInputValue_c;...}
% so that the fitst column will be of the variable names and the second of
% the default values you'd like to assign to them.
%
% vararginInput should contain inputs in the form of a cell array as such:
% {'inputName1',inputValue1,'inputName2',inputValue2,...}
% this is the varargin of the mother function

p = inputParser;
for i = 1:size(defaultInputs,1)
    addParameter(p,defaultInputs{i,1},defaultInputs{i,2});
end
parse(p,vararginInput{:});
Param = p.Results;  % structure containing input parameters of this function
end