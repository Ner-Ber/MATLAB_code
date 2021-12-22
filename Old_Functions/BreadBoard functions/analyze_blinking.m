%% parameters 
bitDepth = 12;
plot_division = 0;
plot_difference = 1;

%% aquire
%---choose folders and files to compare
folders_cell = {};
folder_name = true;
while folder_name
    folder_name = uigetdir('C:\Users\owner\Documents\Shahar_Neri\friction_videos\blinking check\LIGHT_BULB',...
        'CLICK CANCEL WHEN DONE');
    folders_cell = cat(1, folders_cell, folder_name);
end
folders_cell = folders_cell(1:end-1);   %delete last empty cell

%---scan *.tif files in folders selected
image_names = cell(size(folders_cell));
for takeIdx = 1:length(folders_cell)
    listing = dir([folders_cell{takeIdx},'\*.tif']);
    image_names{takeIdx} = {listing.name}';
end

%---create all_frames for movies and create division frames
all_frames_cell = cell(size(image_names));
odd_frames_cell = cell(size(image_names));
even_frames_cell = cell(size(image_names));
for takeIdx = 1:length(image_names)
    all_frames_current  = [];
    for FramIdx = 1:length(image_names{takeIdx})
        currentImage = double(imread(fullfile(...
            folders_cell{takeIdx},image_names{takeIdx}{FramIdx})))./2^(bitDepth);
        all_frames_current = cat(3, all_frames_current, currentImage);
    end
    %---create division frames:
    all_frames_cell{takeIdx} = all_frames_current;
    odd_frames_cell{takeIdx} = all_frames_current(:,:,1:2:end);
    even_frames_cell{takeIdx} = all_frames_current(:,:,2:2:end);    
end

%% analyze division
%---create cells
odd_mean_cell = cellfun(@(a) mean(a,3), odd_frames_cell,'UniformOutput',0);
even_mean_cell = cellfun(@(a) mean(a,3), even_frames_cell,'UniformOutput',0);
division_cell = cellfun(@(a,b) a./b, odd_mean_cell,even_mean_cell,'UniformOutput',0);
%---plot
num_takes = length(division_cell);
if plot_division
    figure;
    for takeIdx = 1:num_takes
        subplot(num_takes,1,takeIdx);
        imagesc(division_cell{takeIdx});
    end
end

%% analyze difference

%---create cells
difference_cell = cell(size(all_frames_cell));
mean_difference_cell = cell(size(all_frames_cell));
for takeIdx = 1:num_takes
    num_odd_frames = size(odd_frames_cell{takeIdx},3);
    num_even_frames = size(even_frames_cell{takeIdx},3);
    diff_in_frames = num_odd_frames - num_even_frames;
    
    if num_odd_frames>num_even_frames
        difference_cell{takeIdx} = ...
            odd_frames_cell{takeIdx}(:,:,1:end-diff_in_frames)...
            - even_frames_cell{takeIdx}(:,:,1:end);
        
    elseif num_odd_frames<num_even_frames
        difference_cell{takeIdx} =...
            odd_frames_cell{takeIdx}(:,:,1:end)...
            - even_frames_cell{takeIdx}(:,:,1:end-diff_in_frames);
        
    else
        difference_cell{takeIdx} =...
            odd_frames_cell{takeIdx}(:,:,:)...
            - even_frames_cell{takeIdx}(:,:,:);
    end
    
    mean_difference_cell{takeIdx} = mean(difference_cell{takeIdx},3);
end


%---find minimal and maximal difference values
all_differences = cell2mat(mean_difference_cell);
min_diff = min(all_differences(:));
max_diff = max(all_differences(:));

%---pick random pixels
% pixels_in_frame = numel(difference_cell{1});
% pix_amount_2analyze = round(0.1*pixels_in_frame);
% pixels_2analyze = randperm(pixels_in_frame, pix_amount_2analyze);


if plot_difference
    frame_means = figure;
    difference_hist = figure;
    for takeIdx = 1:num_takes
        figure(frame_means);
        subplot(num_takes,1,takeIdx);
        imagesc(mean_difference_cell{takeIdx});
        caxis([min_diff max_diff]);
        colorbar;
        
        figure(difference_hist);
        subplot(num_takes,1,takeIdx);
        hist(mean_difference_cell{takeIdx}(:));
    end
end


