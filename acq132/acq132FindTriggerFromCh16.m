function tTrigger=acq132FindTriggerFromCh16(exp_dir)

path=[exp_dir '\acq132_094\stream'];
[ch_094 t]=acq132_stream_read(path);

tTrigger_logical1=logical(([1*(diff(ch_094(:,16))<-2)' 0]));
tTrigger_logical2=logical([0 1*(diff(ch_094(:,16))<-2)']);

tTrigger=(t(tTrigger_logical1)+t(tTrigger_logical2))/2;


