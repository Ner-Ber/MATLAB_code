function blockEdgesCorrdinates = IDT_findEdgePixlesAllEvents(exp_dir, varargin)
    % blockEdgesCorrdinates = IDT_findEdgePixlesAllEvents(exp_dir, varargin)
    %   optional: events,pre_time,intervalTime,post_time,timeBase,smt,lineNum
    %
    %IDT_findEdgePixlesAllEvents will find the events of the block by
    %taking the two outmost results of all events in a current experiment
    
    %% set defaults
    [events,pre_time,intervalTime,post_time,timeBase,smt,lineNum] = setDefaults4function(varargin,...
        'all',0.5e-3,1e-5,4e-3,1,0,'all');
    
    %% set events
    outdata = phantomReadMeta([exp_dir,'\Ph']);
    NumEvents = outdata.NumEvents;
    if strcmpi(events,'all')
        eventsVec = 1:NumEvents;
    else
        eventsVec = events(events<NumEvents & events>1);
    end
    
    %% read row over time
    leftBound = [];
    rightBound = [];
    for i = eventsVec
        try
            ims = phantomGetLines_BigPic(exp_dir,i,-pre_time,intervalTime,post_time,timeBase,smt,'all',lineNum);
        catch
            continue
        end
        RowOverTime_i = ims.lines;
        blockEdgesCorrdinates_i = IDT_findEdgePixlesOfBlock(RowOverTime_i);
        leftBound = cat(1,leftBound,min(blockEdgesCorrdinates_i));
        rightBound = cat(1,rightBound,max(blockEdgesCorrdinates_i));
    end
    
    blockEdgesCorrdinates = [min(leftBound), max(rightBound)];