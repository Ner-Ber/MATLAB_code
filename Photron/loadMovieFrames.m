function [ImageArray, MovieData] = loadMovieFrames(folderPath, frameNumbers)

%--- length of maximal index from frameNumbers
maxImageIndxLength = max(max(ceil(log10(abs(frameNumbers))),1));

MovieData = [];
%--- create list of relevant images file names
ImageArray = [];
folderInfo = dir(folderPath);
for i = 3:length(folderInfo)
    curentFileName = folderInfo(i).name;
    curentFileNum = str2double(curentFileName(end-(4+maxImageIndxLength-1):end-4));
    
    if strcmpi(curentFileName(end-2:end),'cih')
        chiPath = fullfile(folderPath,curentFileName);
        MovieData = Read_cih_file(chiPath(1:end-4));
        
    elseif ismember(curentFileNum,frameNumbers)
        ImportedIm = imread(fullfile(folderPath,curentFileName));
        ImageArray = cat(3,ImageArray,ImportedIm);
    end
end

ImageArray = im2double(ImageArray);