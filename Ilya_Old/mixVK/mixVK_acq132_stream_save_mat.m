function mixVK_acq132_stream_save_mat(exp_dir,file_name,file_name_save)
%exp_dir  '10-48-10' or '10-48-10\cal'
%file_name 'stream'
%some times during the cal needed to save few files. use file_name as
%stream,stream2 and so.

if (nargin<3)
file_name_save=file_name;
end

from=1; %usually you want to save all the data

current_dir=pwd;
%--------------read the data and correct offset units
path=[current_dir '\' exp_dir '\acq132_093\' file_name];
[ch_093,t_093]=acq132_stream_read(path);
path=[current_dir '\' exp_dir '\acq132_094\' file_name];
[ch_094,t_094]=acq132_stream_read(path);

if length(t_093)==length(t_094)
    to=length(t_093);
else if length(t_093)<length(t_094)
        sprintf('the cards arent the same length')
        to=length(t_093);
    else
        sprintf('the cards arent the same length')
        to=length(t_094);
    end
end


%use the commented script to cut the edges correctly, notice that 'end' will make troubles
% if nargin>4
%     from=from-(smt-1)/2;
%     to=to+(smt-1)/2;
% end

ch_093=ch_093(from:to,:);
ch_094=ch_094(from:to,:);
t=t_093(from:to)';

% if length(ch_093(:,1))<length(ch_094(:,1)) %Sometimes the length of the channels is not the same.
% ch_094=ch_094(1:length(ch_093(:,1)),:);
% t=t(1:length(ch_093(:,1)));
% else
% ch_093=ch_093(1:length(ch_094(:,1)),:);
% t=t(1:length(ch_094(:,1)));
% end

acq132_stream=mixVK_acq132_convert_raw_to_data(exp_dir,ch_093,ch_094);
acq132_stream.t=t;
save_path=[current_dir '\' exp_dir '\' file_name_save '.mat'];
save(save_path,'-struct','acq132_stream')

