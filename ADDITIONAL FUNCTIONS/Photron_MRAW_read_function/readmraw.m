function [ImageData] = readmraw(filename,numimgs)
%   readmraw.m 
%   READMRAW Read 16 bit monchrome Photron image data into a structure
%   
%   I=READMRAW('c:\Photron\Filename',[a,b]) loads images a through b from
%   'c:\Photron\Filename' into structure I.  
%
%   Remarks
%   -------
%   This function must be handed the common *.cih and *.mraw file name 
%   and the range of images to be loaded (Matlab can not handle the entire
%   image range for large files).  A file extension is not allowed.
%   This function is intended for monochrome 16 bit *.mraw files only.
%   NOTE: Both the *.cih file and the *.mraw file are utilized
%   Autor: SEP                                Creation Date: June 20,2013
%
%   Examples
%   --------
%   % Load all images
%   I=readmraw('c:\Photron\Moviefile',[0]); 
%   % Load images 10 through 50
%   I=readmraw('c:\Photron\Moviefile',[10,50]);
%   % Load image 10
%   I=readmraw('c:\Photron\Moviefile',[10]);
fid1=fopen(sprintf('%s.cih',filename),'r');
fid2=fopen(sprintf('%s.mraw',filename),'r');
if fid1<1 || fid2<1 display([filename ' filenames could not be found']); ImageData=0; return; end

% Read Header Information
Header=textscan(fid1,'%s','delimiter',':');
Total_Frames=str2double(cell2mat(Header{1}(32)));
if isnan(Total_Frames)==1; 
    Total_Frames=str2double(cell2mat(Header{1}(31)));
    Width=str2double(cell2mat(Header{1}(39)));
    Height=str2double(cell2mat(Header{1}(41)));
    CameraSetup.FrameRate=str2double(cell2mat(Header{1}(23)));
    %Other data from header if desired
else %User Defined Camera Name defined,  shift 1 bit
Width=str2double(cell2mat(Header{1}(40)));
Height=str2double(cell2mat(Header{1}(42)));
CameraSetup.FrameRate=str2double(cell2mat(Header{1}(24)));
%Other data from header if desired
end
Pixels=Width*Height;
fclose(fid1);

% Define Image Range
if numimgs==0               % load all the images 
    first_frame=1;
    frames=Total_Frames;
elseif (length(numimgs)==1) % load a single image
    first_frame=numimgs;
    frames=1;
else                        % load a specified range of images
    first_frame=numimgs(1,1);
    last_frame=numimgs(1,2);
    frames=last_frame-first_frame+1;
end

%Load Images
fseek(fid2,(first_frame-1)*Pixels*2,'bof');
I=zeros(Pixels,frames,'uint16');
for n=1:1:frames
    I(:,n)=(fread(fid2,Pixels,'uint16'));
end
fclose(fid2);
N = [Width Height frames];
ImageData.Images.RawImages=permute(reshape(I,N),[2 1 3]);
ImageData.CameraSetup=CameraSetup;


