%% import images
folderPath = 'C:\Users\owner\Documents\Shahar_Neri\Experiments\2017-11-07\20170711_164056';
frameNumbers = [50:100];
[Images, MovieData] = loadMovieFrames(folderPath, frameNumbers);

%--- crop images:
Images = Images(89:end,:,:);
%--- rotate images
deg2rotation = -5;
figure; imagesc(imrotate(Images(:,:,1),deg2rotation));

rotatedImages = imrotate(Images,deg2rotation);
Mrot = ~imrotate(true(size(Images)),deg2rotation);
rotatedImages(Mrot&~imclearborder(Mrot)) = NaN;

% rotatedImages = imrotate(Images,deg2rotation);

%% differentiate images
ImagesMatSizes = size(rotatedImages);
%     ImagesMatSizes(3) = ImagesMatSizes(3)-1;
%         DiffImages = zeros(ImagesMatSizes);
%
%         for i = 1:ImagesMatSizes(3)
%             DiffImages(:,:,i) = Images(:,:,i+1)-Images(:,:,i);
%         end
%         DiffImages = diff(Images,1,3);

DiffImages = rotatedImages - repmat(rotatedImages(:,:,1),[1,1,size(rotatedImages,3)]);

%     DiffImages = Images;

%% calc differences
% colormap jet
if 0
    hp = {};
    CAxis = [min(DiffImages(:)), max(DiffImages(:))];
    figure;
    for i = 1:ImagesMatSizes(3)
        subplot(ImagesMatSizes(3),1,i)
        imagesc(DiffImages(:,:,i));
        caxis(CAxis);
        hp{i} = get(subplot(ImagesMatSizes(3),1,i),'Position');
    end
    colorbar('Position', [hp{end}(1)+hp{end}(3)+0.03  hp{end}(2)  0.01  (hp{1}(2)+hp{1}(4))-hp{end}(2)]);
    caxis(CAxis);
end

%% create slice of image overtime (heat plot)
if 0
    rowSliceNum = 132;
    ImagesSlice = permute(DiffImages(rowSliceNum,:,:),[3 2 1]);
    
    timeScale = (frameNumbers+MovieData.StartFrame);%*MovieData.ShutterSpeed_s;
    
    figure;
    imagesc(ImagesSlice,'YData',[timeScale(1) timeScale(end)]);
end

%% mean rotated-diff images and heat map
if 1
    meanedImages = mean(DiffImages,1,'omitnan');
    meanedImages = permute(meanedImages, [3 2 1]);
    ColorsMat = MyVaryColor(size(meanedImages,1));
    
    %--- plot row intensity by frame
    figure;
    hold on;
    for i = 1:size(meanedImages,1)
        plot(meanedImages(i,:),'.','Color',ColorsMat(i,:))
    end
    hold off;
    Cbar = colorbar; ylabel(Cbar,'frame num');
    title('meaned row by frame');
    
    %--- plot heat map of differences
    figure; imagesc(meanedImages);
    
end

%% plot rows as  power spectrum
if 0
    Images = Images(5:9,:,:);
    Images = mean(Images,1);
    ColorMat = MyVaryColor(size(Images,3));
    %     row_num = 6;
    figure;
    for rowNumIndx = 1:size(Images,1)
        subplot(size(Images,1),1,rowNumIndx);
        hold all;
        for i = 1:size(Images,3)
            plot(Images(rowNumIndx,:,i),'Color',ColorMat(i,:));
        end
        hold off;
    end
    
end