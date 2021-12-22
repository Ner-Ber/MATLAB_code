function RowOverTimeCell = Movie_readImages(experNum, eventNum)
% create a RowOverTimeCell from the videos.

%% define defaults
exper=my_dir;
disp('');

%% Create the Images matrix, [rows * spatial * time]
disp('Reading images...');
pre_frames=1000;
post_frames=4000;
% [Images,~] = IDT_ReadImages(exper{experNum}, eventNum,pre_frames,1,post_frames,0);      % IDT
[Images,~] = photronReadImages(exper{experNum},eventNum,pre_frames,1,post_frames,0);   %Photron

flip_images =0;
if flip_images ==1
    Images = -Images;
    Images = Images+min(Images(:));
end

%% Make The Row-Over-Time Pictures, each cell for a different row
disp('creating RowOverTimeCell...');
Rows=1:size(Images,1);                 %% choose rows to look at
Mean_RowOverTime=mean(Images,1);
Mean_RowOverTime=permute(Mean_RowOverTime,[3 2 1]);

RowOverTimeCell = {};
for rowIdx = 1:length(Rows)
    RowOverTime = Images(rowIdx,:,:);
    RowOverTime = permute(RowOverTime,[3 2 1]);
    RowOverTimeCell{rowIdx} = RowOverTime;
end

RowOverTimeCell{length(Rows)+1}=Mean_RowOverTime;


end