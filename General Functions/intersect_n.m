function [intersections] = intersect_n(varargin)

if length(varargin)==1 && iscell(varargin)
    VAR = varargin{1};
else
    VAR = varargin;
end

intersections = intersect(VAR{1},VAR{2});
for i=3:length(VAR)
    intersections = intersect(intersections,VAR{i});
end

