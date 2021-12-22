exper = my_dir;
ExpCell = {};
for Exp = 4:6
    EventsCell = {};
    for Ev = 1:30
        try
            CamMetaStruct = CameraMetaAllCams(exper{Exp});
            RowOverTime = phantom_getRowOverTime(exper{Exp}, Ev, CamMetaStruct, 3e-3, 2e-3, 'all');
            EventsCell{Ev} = RowOverTime;
        catch
        end
    end
    ExpCell{Exp} = EventsCell;
end