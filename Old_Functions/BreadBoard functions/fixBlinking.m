function [fixed_frames, difference_mean] = fixBlinking(all_frames, FirstFrames, intensityThreshold)
%fixBlinking will create new frames with averaged-out blinking, according
%to the movie inserted.
% algorithm:
% 		a. Pick frames in beginning of movie where frame doesn't change
% 		b. Pick pixels that are lit above a certain threshold (yet to be decided)
% 		c. Find out which set of frame (odd or even) are in higher intensity (because of blinking)
% 		d. Find mean of difference of the pixels frame stages a & b
% 		e. Subtract the value from d from all frame-set found in c.
% 		f. Export new movie.
% all_frames - a matrix containing frames from a blinking movie. if
% grayscale movie, 3rd dimension of should represent time progress. if RGB
% movie, 3rd dimension should represent RGB layers, 4th- time progress. in
% this case the function will convert frames to YIQ and change intensity
% only. 
% FirstFrames - is a scalar indicating the number of frames (counting from
% the first one) where the frame is still and the difference can be
% calculated.
% fixed_frames - fixed movie with fixed blinking
% difference_mean - is the difference subtracted from half of the frames.
% intensityThreshold - threshold of "lit" pixels, should be <1, >0
% fixBlinking is limited to a blinking between even and odd frames only. 
% the function

%% get frames to work with
    Dimensions = ndims(all_frames);
    if Dimensions==4
        [m n c t] = size(all_frames);
        all_frames_reshaped =  reshape(all_frames, m*t, n, c);
        all_frames_YIQ = reshape(rgb2ntsc(all_frames_reshaped), [m n c t]);
        intensity_frames = permute(all_frames_YIQ(:,:,1,:),[1 2 4 3]);
    else
        intensity_frames = all_frames;
    end
    reference_frames = intensity_frames(:,:,1:FirstFrames);

%% find lit pixels
    litframesThreshold = 0.9;
    %#scalar in range [0,1], indicates fraction of total frames that are required so pixel would count as "lit"

    IntensityThresholdLogical = reference_frames>=intensityThreshold;
    LitPixelsLogical = (sum(IntensityThresholdLogical,3)./FirstFrames)>=litframesThreshold;
    %#a logical frame indicating as pixels that count as lit



    %---devide into odd and even frames
    odd_frames = reference_frames(:,:,1:2:end);
    even_frames = reference_frames(:,:,2:2:end);

%% find higlighted frames and fix
    num_odd_frames = size(odd_frames,3);
    num_even_frames = size(even_frames,3);
    diff_in_frames = num_odd_frames - num_even_frames;

    if num_odd_frames>num_even_frames       %when more odd frames
        difference_frames = ...
            odd_frames(:,:,1:end-diff_in_frames)...
            - even_frames(:,:,1:end);

    elseif num_odd_frames<num_even_frames	%when more even frames
        difference_frames =...
            odd_frames(:,:,1:end)...
            - even_frames(:,:,1:end-diff_in_frames);

    else                                %when same mount of odd and even
        difference_frames =...
            odd_frames(:,:,:)...
            - even_frames(:,:,:);
    end

    difference_mean_frame = mean(difference_frames,3);  %mean difference of each pixel by time
    difference_mean_frame(~LitPixelsLogical) = 0;       %don't change dark pixels
%     difference_mean = mean(difference_mean_frame(LitPixelsLogical));
    
%     difference_mean = mean(difference_frames(find(LitPixelsLogical)));

    fixed_intensity = intensity_frames;
%     fixed_intensity(:,:,1:2:end) = fixed_intensity(:,:,1:2:end) - difference_mean;
    
    fixed_intensity(:,:,1:2:end) = fixed_intensity(:,:,1:2:end)...
        - repmat(difference_mean_frame, [1 1 size(fixed_intensity(:,:,1:2:end),3)]);
    difference_mean = difference_mean_frame;

    if Dimensions==4
        fixed_intensity = permute(fixed_intensity,[1 2 4 3]);
        all_frames_YIQ(:,:,1,:) = fixed_intensity;
        all_frames_YIQ_reshaped =  reshape(all_frames_YIQ, m*t, n, c);
        fixed_frames = reshape(ntsc2rgb(all_frames_YIQ_reshaped), [m n c t]);
    else
        fixed_frames = fixed_intensity;
    end


end