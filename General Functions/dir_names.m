function content_names = dir_names(path)
%content_names = dir_names(path)
%
%dir_names gives a cell array containing the names of files in the 'path'
%given as input

if (nargin==0)
    path=pwd;
end

Content = dir(path);
content_names = {};
for i = 3:size(Content,1)
    content_names{i-2} = Content(i).name;
end


