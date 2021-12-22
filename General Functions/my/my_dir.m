function directoryList=my_dir(path,isdir)
if (nargin==0)
    path=pwd;
    isdir=1;
end
if (nargin==1)
    isdir=1;
end

tmpDirectoryList=dir(path);
directoryListIndex=1;

for j=1:length(tmpDirectoryList)
    if (tmpDirectoryList(j).isdir==isdir && ~strcmp(tmpDirectoryList(j).name,'.') && ~strcmp(tmpDirectoryList(j).name,'..'))
       directoryList{directoryListIndex}=tmpDirectoryList(j).name; 
       directoryListIndex=directoryListIndex+1;
    end
end