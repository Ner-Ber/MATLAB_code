%%%% plot location of phedis systematically
% close all;
% clear;
exper = my_dir;


moveOn = 1;
failNum=0;
ExperNum = 4;
EvNum=17;
PhediStruct={};
ROT_figs = [];
loc_figs = [];
while moveOn
    try
        disp(['*** current event: ',num2str(EvNum),'***'])
        %--- get phedi data
        PhediStruct{EvNum} = Movie_phedi_from_folder_2_data_defByName(ExperNum,EvNum,...
            'preRowsTime',3e-3,'postRowsTime',3e-3,'prePhediTime',2e-3,'postPhediTime',2e-3,'QuickAndDirty',1);
        
        %--- plot the big pic ROT
        ROT_figs(EvNum) = figure;
        ax(1) = subplot(2,1,1);
        IDT_PlotRowOverTime(PhediStruct{EvNum}.BigPicRotStruct);
        %--- plot the phedi ROT
        ax(2) = subplot(2,1,2);
        axis off;
        
        Movie_phedi_plotWithSg(PhediStruct{EvNum},1,'sep',{'ROT'});
        ROT_phedi_handle = gcf;
        h = get(ROT_phedi_handle,'Children');
        newh = copyobj(h,ROT_figs(EvNum));
        for j=1:length(newh)
            posnewh = get(newh(j),'Position');
            possub  = get(ax(1),'Position');
            set(newh(j),'Position',...
                [posnewh(1) posnewh(2) posnewh(3) possub(4)])
        end
        close(ROT_phedi_handle);
        
        Movie_phedi_plotWithSg(PhediStruct{EvNum},1,'sep',{'loc'});
        loc_figs(EvNum) = gcf;
        
        
        EvNum = EvNum+1;
    catch
        EvNum = EvNum+1;
        failNum=failNum+1;
        if failNum>=5
            break
        end
    end
    
end