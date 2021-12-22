function [CAM_Meta] = CameraMetaAllCams(ExperimentPath,FolderName)

folderInExperDir = my_dir(ExperimentPath);

if nargin==1
    if sum(strcmpi('Ph',folderInExperDir))>0
        CAM_Meta = phantomReadMeta([ExperimentPath,'\Ph']); %Phantom
        CAM_Meta.CameraName = 'Ph_small';
    elseif sum(strcmpi('IDT',folderInExperDir))>0
        CAM_Meta = IDT_ReadMeta([ExperimentPath,'\IDT']);   %IDT
        CAM_Meta.CameraName = 'IDT';
    elseif sum(strcmpi('photron',folderInExperDir))>0
        CAM_Meta = photronReadMeta([ExperimentPath,'\photron']);   %Photron
        CAM_Meta.CameraName = 'Photron';
    end
else
    switch FolderName
        case 'Ph'
            CAM_Meta = phantomReadMeta([ExperimentPath,'\Ph']); %Phantom
            CAM_Meta.CameraName = 'Ph_small';
        case 'PhBig'
            CAM_Meta = phantomReadMeta([ExperimentPath,'\PhBig']); %Phantom
            CAM_Meta.CameraName = 'Ph_big';
        case 'IDT'
            CAM_Meta = IDT_ReadMeta([ExperimentPath,'\IDT']);   %IDT
            CAM_Meta.CameraName = 'IDT';
        case 'photron'
            CAM_Meta = photronReadMeta([ExperimentPath,'\photron']);   %Photron
            CAM_Meta.CameraName = 'Photron';
    end
    
end

end