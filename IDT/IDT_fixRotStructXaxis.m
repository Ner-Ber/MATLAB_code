function outputROT = IDT_fixRotStructXaxis(RotStruct)

blockLength = 0.15;     % default
[~,LeftBound] = min(abs(RotStruct.x-0));
[~, RightBound] = min(abs(RotStruct.x-blockLength));

%% choose bounds
coor_cell = {};
measureFig = figure;
imagesc(RotStruct.DataMatNorm);
title('mark left or right bound. only last click will count');
caxis([0.5 1.1]);
[~, leftTick] = min(abs(RotStruct.x));
[~, rightTick] = min(abs(RotStruct.x-blockLength));
set(gca, 'XTick', [leftTick rightTick], 'XTickLabel', {num2str(0),num2str(blockLength)});

%--- create buttons for UI
uicontrol('Style', 'pushbutton', 'String', 'left boundary',...
    'Position', [10 20 80 40],...
    'Callback',@MarkLeft);

uicontrol('Style', 'pushbutton', 'String', 'right boundary',...
    'Position', [10 60 80 40],...
    'Callback',@MarkRight);

uicontrol('Style', 'pushbutton', 'String', 'Change block from 0.15m',...
    'Position', [10 100 80 40],...
    'Callback',@ChangeLength);

%--- terminate figure containing UI
uiwait(gcf);         % wait for command and uiresume

%% create new fixed structure
newX = nan(size(RotStruct.x));
newX(round(LeftBound(end))) = 0;
newX(round(RightBound(end))) = blockLength;
newX = interp1(round([LeftBound(end), RightBound(end)]),[0 blockLength],1:length(newX),'linear','extrap');

outputROT = RotStruct;
outputROT.x = newX;

%% nested functions
    function MarkLeft(src,event)
        [LeftBound,~] = ginput;
    end

    function MarkRight(src,event)
        [RightBound,~] = ginput;
    end

    function ChangeLength(src,event)
        prompt = 'Enter block length:';
        dlg_title = 'change block length';
        num_lines = 1;
        defaultans = {'0.15'};
        blockLength = inputdlg(prompt,dlg_title,num_lines,defaultans);
        blockLength = str2double(blockLength);
    end
end