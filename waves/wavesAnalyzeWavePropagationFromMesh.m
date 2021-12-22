function eventsCell = wavesAnalyzeWavePropagationFromMesh(DataStructCell)
% analyzePhedi_detrend_and_strech


%% create buttons for UI
done = 0;
MeshFig = figure;
eventsCell = {};
i = 1;
eventsCell{i} = [];
plotImagesc(i);



uicontrol('Style', 'pushbutton', 'String', 'Next event',...
    'Position', [10 20 80 40],...
    'Callback',@NextEv);

uicontrol('Style', 'pushbutton', 'String', 'Previous event',...
    'Position', [10 60 80 40],...
    'Callback',@PrevEv);

uicontrol('Style', 'pushbutton', 'String', 'jump 2 Ev',...
    'Position', [10 100 80 40],...
    'Callback',@jump2Ev);

uicontrol('Style', 'pushbutton', 'String', 'measure line',...
    'Position', [10 140 80 40],...
    'Callback',@measureLine);
uicontrol('Style', 'pushbutton', 'String', 'Done',...
    'Position', [10 180 80 40],...
    'Callback',@Done);


%--- terminate figure containing UI
waitfor(MeshFig);         % wait for command and uiresume

%% nested functions
    function NextEv(src,event)
        i = i+1;
        plotImagesc(i)
        eventsCell{i} = [];
    end

    function PrevEv(src,event)
        i = i-1;
        plotImagesc(i)
        eventsCell{i} = [];
    end

    function jump2Ev(src,event)
        prompt = 'Event number';
        dlg_title = 'Jump to event';
        num_lines = 1;
        defaultans = {'4'};
        Event2Jump = inputdlg(prompt,dlg_title,num_lines,defaultans);
        Event2Jump = str2double(Event2Jump);
        i = Event2Jump;
        plotImagesc(i);
        eventsCell{i} = [];
    end

    function measureLine(src,event)
%         position_cell = {};
        h = imline;
        position = wait(h);
        hold on; line(position(:,1),position(:,2),'LineWidth',1.2,'Color','r');
        eventsCell{i} = cat(3,eventsCell{i},position);
    end

    function Done(src,event)
        close(MeshFig); 
    end

    function plotImagesc(i)
        try
            DataStruct = DataStructCell{i};
            %-- strech phedi data:
            PhediLocation = DataStruct.PhediData.PhediLocation;
            timeVec = DataStruct.PhediData.timeVec;
            PhediLocation_raw_stretch = stretchAmpMinus1to1(PhediLocation);
            slopeLogical = DataStruct.PhediData.slopeIncline>0;
            %--- create data for RowOverTime
            x = [timeVec(1), timeVec(end)];
            y = [1,size(PhediLocation_raw_stretch(:,slopeLogical),2)];
            FN = DataStruct.SgData.N;
            FS = DataStruct.SgData.F;
            ExpState = ['Nmean=',num2str(round(mean(FN))),'N1=',num2str(round(mean(FN(1:100)))),' Nend=',num2str(round(mean(FN(end-100:end)))),...
                '  Fmean=',num2str(round(mean(FS))),'F1=',num2str(round(mean(FS(1:100)))),' Fend=',num2str(round(mean(FS(end-100:end))))];
            %--- plot
            IMGSC = imagesc(x,y,PhediLocation_raw_stretch(:,slopeLogical)');
            colormap(jet);
            ax = gca;
            ax.YDir = 'normal';
            title({'mesh of raw data normalized', DataStruct.BigPicRotStruct.details,ExpState});
            xlabel('t [s]');
            ylabel('position on block [arb]');
        catch
            plot(nan,nan);
            axis([-1 1 -1 1]);
            text(0,0,['could not plot imagesc for event ',num2str(i)]);
            warning(['could not plot imagesc for event ',num2str(i)]);
        end
        
    end
end