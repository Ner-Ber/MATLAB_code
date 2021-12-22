function ImageData = get_cih_data(folderPath)

switch folderPath
    case ''
        [file_name,folderPath,~] = uigetfile({'*.cih','cih file';...
        '*.*','All Files' },'select images to find edges (**same batch**)','MultiSelect','on');
    otherwise
        listing = dir([folderPath,'\*.cih']);
        file_name = {listing.name}';
        file_name = file_name{:};
end
file_name = file_name(1:end-4);

ImageData = Read_cih_file(fullfile(folderPath,file_name));
end