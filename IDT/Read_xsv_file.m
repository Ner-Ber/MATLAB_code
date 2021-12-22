function [xsvData] = Read_xsv_file(filename)


%---specify data to aquire
wantedFieldNames = {'Date', 'Time', 'FrameRate', 'FrameT',...
    'EffectiveBitDepth','ImageWidth','ImageHeight', 'NumIms',...
    'StartFrame'};

%% open file to read
fid1=fopen(filename,'r');
if fid1<1
    display([filename ' filenames could not be found']);
    xsvData=0;
    return;
end

%% organize data
%---Read Header Information
Header=textscan(fid1,'%s','delimiter','=');
Header = Header{1};
CategoryIdx = regexp(Header,'\[\w*[a-z]\w*\]');
CategoryIdx = ~cellfun(@isempty,CategoryIdx);
if nnz(CategoryIdx)>1
    CategoryIdxLoc = find(CategoryIdx);
    Header = Header(2:(CategoryIdxLoc(2)-1));
else
    Header = Header(2:(CategoryIdxLoc(2)-1));
end

Names = Header(1:2:end);
Values = Header(2:2:end);

Names = Names(1:min(length(Values),length(Names)));
Values = Values(1:min(length(Values),length(Names)));

%--- diminish empty spaces
Names = strrep(Names, ' ', '');
Names = strrep(Names, '-', '');

%% turn numbers into doubles:
NumericValues = Values;
for i = 1:length(Values)
    NumericValues{i} = str2double(Values{i});
end
NumericValues(cellfun(@isnan,NumericValues)) = Values(cellfun(@isnan,NumericValues));

%% write structure
xsvData = cell2struct(NumericValues',Names,2);

%% add wanted fields
Data2Add = {[num2str(xsvData.Day),'/',num2str(xsvData.Month),'/',num2str(xsvData.Year)],...
    [num2str(xsvData.Hour),':',num2str(xsvData.Minute)],...
    xsvData.Rate, xsvData.ExposureNs,xsvData.PixelDepth,...
    xsvData.ROIW,xsvData.ROIH, xsvData.BROCLenght,...
    xsvData.PerTrigFrms};

for i = 1:length(wantedFieldNames)
    xsvData.(wantedFieldNames{i}) = Data2Add{i};
end


%% done 
fclose(fid1);

