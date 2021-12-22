function [coors_over_frames] = Movie_follow_point_in_movie(Images)

coors_over_frames = {};
figureFrame = figure; %('PropertyName','figureFrame');

%% CREATE PUSHBUTTONS
%--- done marking this trajectory button:
DoneTrajButton = uicontrol('Style', 'pushbutton', 'String', 'New Trajectory',...
    'Position', [20 20 100 40],...
    'Callback', @DoneTrajectoryMarking);
%--- finished marking this movie
DoneMovieButton = uicontrol('Style', 'pushbutton', 'String', 'Done Marking',...
    'Position', [120 20 100 40],...
    'Callback', @DoneMovieMarking);
%--- Jump to frame
Jump2FrameButton = uicontrol('Style', 'pushbutton', 'String', 'Jump 2 Frame',...
    'Position', [400 20 100 40],...
    'Callback', @Jump2Frame);
%--- Perv frame
PrevFrameButton = uicontrol('Style', 'pushbutton', 'String', 'Prev Frame',...
    'Position', [300 20 100 40],...
    'Callback', @PrvFrame);
%--- Nxt frame
NxtFrameButton = uicontrol('Style', 'pushbutton', 'String', 'Nxt Frame',...
    'Position', [500 20 100 40],...
    'Callback', @NxtFrame);
%--- reset zoom
RetainZoomButton = uicontrol('Style', 'pushbutton', 'String', 'Retain Zoom',...
    'Position', [700 20 100 40],...
    'Callback', @RetainZoom);
%--- frames to jump
frames2JumpButton = uicontrol('Style', 'pushbutton', 'String', 'Frames 2 Jump',...
    'Position', [1000 20 100 40],...
    'Callback', @Frames2Jump);



%% gather data:
trajectory_idx = 1;
FramesJump = 1;
ChangeFrameFlag = 0;
JumpFrameFlag = 0;
NxtFrameFlag = 0;
PrvFrameFlag = 0;
DoneTrajFlag = 0;
DoneMovieFlag = 0;
RetainZoomFlag = 0;

while ~DoneMovieFlag
    coors_over_frames{trajectory_idx} = [];
    
    figureHandle = figure(figureFrame);
    i = 1;
    while i <= size(Images,3)
        %--- display image and measurment details
        imshow(Images(:,:,i));
        title({['Trajectory Length:  ',num2str(size(coors_over_frames{trajectory_idx},1))],...
            ['current Frame:  ',num2str(i)],...
            ['Trajectory Num:   ', num2str(trajectory_idx)]});
        set(gca,'Position',[0 0 1 1])
        
        %--- display in wanted zoom:
        if i==1
            ax = figureHandle.CurrentAxes;
            AxesSetup = [ax.XLim ax.YLim];
        elseif RetainZoomFlag
            axis image
            set(gca,'Position',[0 0 1 1])
            RetainZoomFlag = 0;
        else
            axis(AxesSetup);
        end
        
        %--- display current marking of point:
        if size(coors_over_frames{trajectory_idx},1)>=i
            hold on;
            plot(coors_over_frames{trajectory_idx}(i,1),...
                coors_over_frames{trajectory_idx}(i,2),...
                'rx');
        end
        
        %--- make a measurment:
        h = impoint(gca);
        p = wait(h);
        
        %--- get current zoom:
        ax = figureHandle.CurrentAxes;
        AxesSetup = [ax.XLim ax.YLim];
        
        %--- break if done with measurment
        if DoneTrajFlag || DoneMovieFlag
            break
        end
        
        %--- continue to retain zoom
        if RetainZoomFlag
            continue
        end
        
        %--- change frame by usr request:
        if ChangeFrameFlag
            %----------
            if JumpFrameFlag
                if frameNum>size(Images,3)
                    InvalidFrameDialog();
                else
                    i = frameNum;
                end
            %----------
            elseif NxtFrameFlag
                i = i+FramesJump;
            %----------
            elseif PrvFrameFlag
                i = i-FramesJump;
            end
            %----------
            ChangeFrameFlag = 0;
            JumpFrameFlag = 0;
            NxtFrameFlag = 0;
            PrvFrameFlag = 0;
            continue
        end
        %--- save measurment:
        coors_over_frames{trajectory_idx}(i,:) = p;
        i = i+FramesJump;
    end
    trajectory_idx = trajectory_idx+1;
    DoneTrajFlag = 0;
end
close(figureFrame);

%% nested functions
    function DoneTrajectoryMarking(src,event)
        DoneTrajFlag = 1;
        return
    end


    function DoneMovieMarking(src,event)
        DoneMovieFlag = 1;
        return
    end


    function Jump2Frame(src,event)
        UsrAnswer = inputdlg('Frame t jump to:','Jump 2 frame');
        frameNum = str2double(UsrAnswer);
        ChangeFrameFlag = 1;
        JumpFrameFlag = 1;
    end

    function PrvFrame(src,event)
        ChangeFrameFlag = 1;
        PrvFrameFlag = 1;
    end

    function NxtFrame(src,event)
        ChangeFrameFlag = 1;
        NxtFrameFlag = 1;
    end

    function RetainZoom(src,event)
        RetainZoomFlag = 1;
    end

    function Frames2Jump(src,event)
        prompt = 'Enter number of frames to jump each sample or click:';
        defaultInput = num2str(FramesJump);
        answer = inputdlg(prompt);
        FramesJump = str2double(answer);
    end

    
    function InvalidFrameDialog()
        d = dialog('Position',[300 300 250 150]);
        
        txt = uicontrol('Parent',d,...
            'Style','text',...
            'Position',[20 80 210 40],...
            'String','Frame number inserted exceeds number of total frames');
        
        btn = uicontrol('Parent',d,...
            'Position',[85 20 70 25],...
            'String','Oops',...
            'Callback','delete(gcf)');
    end


end



