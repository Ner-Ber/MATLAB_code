function [figureHandle, HeatMap, RowOverTime] = plotPixelRowOverTime2(Num_of_lines,First_line,pattern_tracing,video_num)
%plotPixelRowOverTime will plot the intensity change of a certain row in
%the images inserted as function of time.
% Num_of_lines - the number of lines you want to average and annalize.   
% First_line - the row number you want to plot, where row1 is top row of
% frame. if this are more than one line, then enter the first

% In case you want to annalize the 'pattern_tracing' movies:
% put 1 in the "pattern_tracing"
% choose video_num - from 1 to 10
%
% starting_frame - the first frame you want to show. the events happen at
% ~12700
% num_of_frame - the number you want to plot. In this case, every 256
% frames is a msec

fps=256000;
pixel_size=3e-6;
starting_frame=12500;
num_of_frame=fps/100+1000;

%% ---load images

switch pattern_tracing
    case 0
        [image_names,figures_folder_path,~] = uigetfile({'*.jpg;*.tif;*.png;*.gif;*.ibw','All Image Files';...
        '*.*','All Files' },'select images to find edges (**same batch**)','MultiSelect','on');
    case 1
        fmt='.tif';
        main_folder='C:\Users\owner\Documents\Shahar_Neri\Pattern_Tracing\sampling';
        figures_folder_path= [main_folder '\sample_' num2str(video_num) '\PatternTrace' num2str(video_num) '_C001H001S0001'];

        image_names={};
        for i=1:num_of_frame
            image_names{i}=['PatternTrace' num2str(video_num) '_C001H001S00010' num2str(starting_frame-1+i) fmt];
        end

%% ---create matrix to plot (containing all relevant rows)
RowOverTime = [];
for FramIdx = 1:length(image_names)
    currentFrame = imread(fullfile(figures_folder_path,image_names{FramIdx}));
    
    all_lines=0;
    for i=First_line:First_line+Num_of_lines-1
      all_lines=all_lines+double(currentFrame(i,:));
    end
    Av_line=all_lines./Num_of_lines;
    
    RowOverTime = cat(1, RowOverTime,  Av_line);
    
end


%% ---Normalize Image
% RowOverTime = RowOverTime./max(RowOverTime(:)); % normalized by maximun intensity
%RowOverTime = RowOverTime - repmat(RowOverTime(1,:), size(RowOverTime,1), 1); % normalized by first row
% RowOverTime(isnan(RowOverTime)) = 1; %NaN are probably non-changing pixels
% RowOverTime(isinf(RowOverTime)) = max(max(RowOverTime(RowOverTime~=inf)));

%% save the matrix

save([main_folder '\sample_' num2str(video_num) '\RowOverTime' num2str(video_num) '.mat'],'RowOverTime');


%% plot the result

figureHandle = figure;

MatSize = size(RowOverTime);
frameDuration = 1/fps;
factor=1e3;
x = linspace(0,MatSize(2)*pixel_size,MatSize(2))*factor;
y = linspace(0,MatSize(1)*frameDuration,MatSize(1))*factor;
HeatMap = imagesc('XData',x,'YData',y,'CData',RowOverTime);
xlim([0 x(end)]);  ylim([0 y(end)]);
xlabel('X [mm]');
ylabel('Time [ms]');

% figureHandle = figure;
% plotHandle = imagesc(RowOverTime);
% colorbar;
% ylim([startFrame FramIdx]);

end