function  [Handles, DiffImages] = photronMovieDifference(Images, FrameNums, difference)
%photronMovieDerivative will create a plot of the differences between
%frames.
%
% Arguments:
% Images - a 3D matrix containing the images from the video. Dimensions 1
%   and 2 are hieight and width, 3rd dimension is frames over time.
% FrameNums - frames indecies (index if the 3 dimension) to consider in
%   matching features.
% difference - a string specifying weather to create difference between
% two consecutive frames of between all frames to the first one. Should be
% specified as 'derivative' of 'first'.
%
% Returns:
% Handles - a structure containig handles of the figure and plots.
% MovieDerivative - a 3D matrix containing differences between the
%   consecutive images from the video


SelectedImages = Images(:,:,FrameNums);
ImagesMatSizes = size(SelectedImages);

if strcmpi(difference, 'derivative')
    ImagesMatSizes(3) = ImagesMatSizes(3)-1;
    DiffImages = zeros(ImagesMatSizes);
    for i = 1:ImagesMatSizes(3)
        DiffImages(:,:,i) = SelectedImages(:,:,i+1)-SelectedImages(:,:,i);
    end
elseif strcmpi(difference, 'first')
    DiffImages = SelectedImages - repmat(SelectedImages(:,:,1),1,1,ImagesMatSizes(3));
end

%% plot differences
CAxis = [min(DiffImages(:)) max(DiffImages(:))];

Handles = struct;
Handles.figureHandle = figure;
Handles.subplots = {};
figure;
for i = 1:ImagesMatSizes(3)
    Handles.subplots{i} = subplot(ImagesMatSizes(3),1,i);
    imagesc(DiffImages(:,:,i));
    caxis(CAxis);
end
colorbar('Position',...
    [Handles.subplots{end}.Position(1)+Handles.subplots{end}.Position(3)+0.03...
    Handles.subplots{end}.Position(2)...
    0.01  ...
    (Handles.subplots{1}.Position(2)+Handles.subplots{1}.Position(4))-Handles.subplots{end}.Position(2)]);
caxis(CAxis);



end