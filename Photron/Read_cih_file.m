function [ImageData] = Read_cih_file(filename)


%---specify data to aquire
Data2Aquire = {'Date', 'Time', 'Record Rate(fps)', 'Shutter Speed(s)',...
    'EffectiveBit Depth','Image Width','Image Height', 'Total Frame',...
    'Start Frame'};
fieldNames = {'Date', 'Time', 'FrameRate', 'FrameT',...
    'EffectiveBitDepth','ImageWidth','ImageHeight', 'NumIms',...
    'StartFrame'};

%---open file to read
fid1=fopen(sprintf('%s.cih',filename),'r');
if fid1<1
    display([filename ' filenames could not be found']);
    ImageData=0;
    return;
end

%---Read Header Information
Header=textscan(fid1,'%s','delimiter',':');
Header = Header{1};

% IndexHeader = strfind(Header, '#');
% Index = find(not(cellfun('isempty', IndexHeader)));
% sectionBegin = Index+1;
% sectionEnd = Index-1;
% sectionRanges = [sectionBegin, [sectionEnd(2:end);length(Header)]];

DataCell = {};
for i = 1:length(Data2Aquire)
    currentData = Data2Aquire{i};
    titlePosition = find(strcmp(strtrim(Header),currentData));
    if ~isempty(strfind(currentData,'Time'))   % data is of time format
        DataCell{i} = [Header{titlePosition+1},':',Header{titlePosition+2}];
    elseif ~isempty(strfind(currentData,'Date'))   % data is of date format
        DataCell{i} = Header{titlePosition+1};
    else                                    % data is of time or date format
        DataCell{i} = str2num(Header{titlePosition+1});
    end
    
end
ImageData = cell2struct(DataCell,fieldNames,2);

fclose(fid1);

