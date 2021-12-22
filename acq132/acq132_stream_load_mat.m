function stream=acq132_stream_load_mat(filePath,Tstart,Tinterval,Tend,smt,varargin)
%stream=acq132_stream_load_mat(exp_dir,from,to,varargin)
%exp_dir is the full path. example:  '19-58-52\cal\stream2.mat'

%the function loads the varaibles indicated by 'varargin' (str), then cuts the data.
%from & to in sec
%varargin - strings of variable names. example: 'Sxy','F' 
%first smoothes the data then subsamples it.

%--load the m file
path=filePath;
varargin=['t' varargin];
[stream]=load(path,varargin{:});

%----cut the data edges
if strcmp(Tstart,'start')
    Tstart=stream.t(1);
end
if strcmp(Tend,'end')
    Tend=stream.t(end);
end
if strcmp(Tinterval,'min')
    intervalInd=1;
else
    intervalInd=floor(Tinterval/(stream.t(2)-stream.t(1)));
end

[~,startInd]=min(abs(Tstart-stream.t));
[~,endInd]=min(abs(Tend-stream.t));
t_length_before_cut=length(stream.t(:,1));

for j=1:length(varargin)
    if length(stream.(varargin{j})(:,1))==t_length_before_cut; %cut the time dependent variables
       stream.(varargin{j})=my_smooth(stream.(varargin{j})(startInd:1:endInd,:),smt);
        %stream.(varargin{j})=stream.(varargin{j})(startInd:intervalInd:endInd,:);
        stream.(varargin{j})=stream.(varargin{j})(1:intervalInd:end,:);
    end
end

% assign date and exp
path=pwd;
lastSlashPosition = find(path == '\', 1, 'last');

stream.Path=path;
stream.Date = path(lastSlashPosition+1:end);
stream.exp=filePath(1:8);

%--- assignin the acq132_out to the caller
if nargout==0
    field_names=fieldnames(stream);
    for j=1:length(field_names)
        assignin('caller',field_names{j},stream.(field_names{j}));
    end
    clear stream;
end


