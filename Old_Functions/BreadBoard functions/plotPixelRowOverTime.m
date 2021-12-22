function [figureHandle, plotHandle, RowOverTime] = plotPixelRowOverTime(rowNum)
%plotPixelRowOverTime will plot the intensity change of a certain row in
%the images inserted as function of time.
% rowNum - row number wanted to plot, where row1 is top row of frame.
% startFrame - first frame number of movie

%---load images
[image_names,figures_folder_path,~] = uigetfile({'*.jpg;*.tif;*.png;*.gif;*.ibw','All Image Files';...
'*.*','All Files' },'select images to find edges (**same batch**)','MultiSelect','on');

%---create matrix to plot (containing all relevant rows)
RowOverTime = [];
for FramIdx = 1:length(image_names)
    currentFrame = imread(fullfile(figures_folder_path,image_names{FramIdx}));
    RowOverTime = cat(1, RowOverTime, double(currentFrame(rowNum,:)));
end




%---Normalize Image
% RowOverTime = RowOverTime./max(RowOverTime(:)); % normalized by maximun intensity
%RowOverTime = RowOverTime - repmat(RowOverTime(1,:), size(RowOverTime,1), 1); % normalized by first row
% RowOverTime(isnan(RowOverTime)) = 1; %NaN are probably non-changing pixels
% RowOverTime(isinf(RowOverTime)) = max(max(RowOverTime(RowOverTime~=inf)));


%---plot the result
figureHandle = figure;
plotHandle = imagesc(RowOverTime);
colorbar;
% ylim([startFrame FramIdx]);



end