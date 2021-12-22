function mixVK_acq132_stream_plot_cont(plot_conf,use_offset)
%plot_conf - plot configuration 'exp' ,'ds', 'sg cal'
%use_offset = 'none' ->no offset correction, otherwise, if
%c:\frics\tmp_cooking\cal exists --> offset will be taken from there,if
%not,from the begining of the stream file.

if nargin<2
    use_offset=0;
end
set(0,'DefaultFigureWindowStyle','docked');

path_Master='C:\Frics\tmp\acq132_093\stream';
path_Slave='C:\Frics\tmp\acq132_094\stream';

num_sampl=3000;
byts_to_plot=num_sampl*2*33;
jump=1; %subsampling the data
average_num=100; %the number of samples for shear and normal mean
Ds_offset=1.01;

exp_dir='c:\Frics\tmp';
%------------Read exp_details.txt
exp_details=mixVK_expDetailsRead(exp_dir);

factor_Master=exp_details.factor_Master;
factor_Slave=exp_details.factor_Slave;
sg_Master_order=exp_details.sg_Master_order;
sg_Slave_order=exp_details.sg_Slave_order;
x_sg=exp_details.x_sg;
%MasterVref=exp_details.MasterVref;
MasterCurrent=exp_details.MasterCurrent;
MasterResistances=exp_details.MasterResistances;
factor_Master([1:15 17:31])=factor_Master([1:15 17:31]).*MasterResistances*MasterCurrent*1e-6*132;
%-----------------offset calculation (no offset for external channels)
cal_dir='c:\Frics\tmp\cal'; % if cal directory exists offset correction taken from there, otherwise from the first event
%plotbrowser('off');
%drawnow
if strcmp(use_offset,'none')
    ch_offset_Master=zeros(1,32);
    ch_offset_Slave=zeros(1,32);
    path_format_Master='C:\Frics\format093';%need to be changed
    path_format_Slave='C:\Frics\format094';
    
else if  exist(cal_dir,'dir')==7
        path='C:\Frics\tmp\cal\acq132_093\multivent';
        ch_Master=acq132_event_read(path,1);%using event1 for offset correction
        path='C:\Frics\tmp\cal\acq132_094\multivent';
        ch_Slave=acq132_event_read(path,1);%using event1 for offset correction
        path_format_Master='C:\Frics\tmp\cal\acq132_093\multivent\01.COOKED\format';
        path_format_Slave='C:\Frics\tmp\cal\acq132_094\multivent\01.COOKED\format';
    else
        %--------use start of 'stream' for mean values calculation
        
        if exist('Z:\01.COOKED\format')
            path_format_Master='Z:\01.COOKED\format';
            path_format_Slave='Y:\01.COOKED\format';
        else
            path_format_Master='C:\Frics\format093';%need to be changed
            path_format_Slave='C:\Frics\format094';
        end
        
        norm_coeff_Master=format_read(path_format_Master);
        norm_coeff_Slave=format_read(path_format_Slave);
        ch_offset_09=zeros(1,32);
        factor_09=ones(1,32);
        
        ch_multiplex=multiplex_read(path_Master,byts_to_plot,'bof');
        ch_Master=multiplex_decoding(ch_multiplex,norm_coeff_Master,ch_offset_09,factor_09,jump);
        ch_multiplex=multiplex_read(path_Slave,byts_to_plot,'bof');
        ch_Slave=multiplex_decoding(ch_multiplex,norm_coeff_Slave,ch_offset_09,factor_09,jump);
        
        %--------use event1 for mean values calculation
        %                 path='z:\';
        %                 ch_093=acq132_event_read(path,1);%using event1 for offset correction
        %                 path='y:\';
        %                 ch_094=acq132_event_read(path,1);%using event1 for offset correction
        %                 path_format_093='Z:\01.COOKED\format';
        %                 path_format_094='Y:\01.COOKED\format';
        %
    end
    ch_offset_Master=mean(ch_Master,1);
    ch_offset_Slave=mean(ch_Slave,1);
    ch_offset_Slave(1,29:32)=0;%no offset for external channels
    clear t ch_Master ch_Slave;
end


norm_coeff_Master=format_read(path_format_Master);
norm_coeff_Slave=format_read(path_format_Slave);

%------------The plot window is constant size

if strcmp(plot_conf,'exp') || strcmp(plot_conf,'ds')
    
    for k=1:3000
        ch_multiplex=multiplex_read(path_Master,byts_to_plot);
        ch_Master=multiplex_decoding(ch_multiplex,norm_coeff_Master,ch_offset_Master,factor_Master,jump);
        ch_multiplex=multiplex_read(path_Slave,byts_to_plot);
        ch_Slave=multiplex_decoding(ch_multiplex,norm_coeff_Slave,ch_offset_Slave,factor_Slave,jump);
        
        
        sg=[mean(ch_Master(1:average_num,sg_Master_order)) mean(ch_Slave(1:average_num,sg_Slave_order))]; %mean ordered strain gages
        [~,Syy,Sxy,~,Uyy,Uxy]=mixVK_calculate_stress_strain(sg,exp_details);
        
        N=ch_Master(:,16);
        F=ch_Master(:,32);
        
        if strcmp(plot_conf,'exp')
            %----------plot stress NF F
            
            subplot(3,2,1);
            plot(N,'.-')
            set(gca,'FontSize',16,'FontWeight','bold')
            title('Fn');
            subplot(3,2,2);
            plot(F,'.-');
            set(gca,'FontSize',16,'FontWeight','bold')
            title('Fs');
            subplot(3,2,3);
            plot(x_sg([1:9 11:end]),-Uyy([1:9 11:end]),'.-');
            set(gca,'FontSize',16,'FontWeight','bold')
            title('Uyy ');
            subplot(3,2,4);
            plot(x_sg,Uxy,'.-');
            set(gca,'FontSize',16,'FontWeight','bold')
            title('Uxy ');
            
            
            %plot(x_sg,(Sxy./Syy),'.-');
            %plot(ch_094(:,17)+Ds_offset,'.-');
            %hold all
            %------plotyy
            %         if k>1
            %             delete(h(1));
            %         end
            %             subplot(3,2,5);
            %             h=plotyy(1:length(ch_094(:,17)),ch_094(:,17)+Ds_offset,1:length(ch_094(:,29)),ch_094(:,29)+Ds_offset);
            
            subplot(3,2,5);
            plot(x_sg,Sxy./Syy,'.-');
            xlim([0 200]);
            ylim([0.4 1]);
            set(gca,'FontSize',16,'FontWeight','bold')
            title('Sxy/Syy ');
            
%             plot(ch_Slave(:,31)+Ds_offset,'.-');
%             set(gca,'FontSize',16,'FontWeight','bold')
%             hold off
%             title('ch-31');
            
            subplot(3,2,6);
            plot(ch_Slave(:,32),'.-');
            hold off
            set(gca,'FontSize',16,'FontWeight','bold')
            title('ch-32');
            
        elseif strcmp(plot_conf,'ds')
            %--------Disp. calibartion Plot
            subplot(2,2,4);
            plot(N,'.-')
            set(gca,'FontSize',16,'FontWeight','bold')
            title('Fn');
            
            %             subplot(3,2,2);
            %             plot(ch_094(:,17)+Ds_offset,'.-');
            %             set(gca,'FontSize',16,'FontWeight','bold')
            %             title('17');
            %
            %             subplot(3,2,3);
            %             plot(ch_094(:,29)+Ds_offset,'.-');
            %             set(gca,'FontSize',16,'FontWeight','bold')
            %             title('29');
            
            subplot(2,2,1);
            plot(ch_Slave(:,32),'.-');
            set(gca,'FontSize',16,'FontWeight','bold')
            axis('tight')
            title('32');
%             subplot(2,2,1);
%             plot(F,'.-')
%             set(gca,'FontSize',16,'FontWeight','bold')
%             title('Fs');
%             
            subplot(2,2,2);
            plot(ch_Slave(:,30)+Ds_offset,'.-');
            set(gca,'FontSize',16,'FontWeight','bold')
            axis('tight')
            title('30');
            
            subplot(2,2,3);
            plot(ch_Slave(:,31)+Ds_offset,'.-');
            set(gca,'FontSize',16,'FontWeight','bold')
            axis ('tight')
            title('31');
            
%             subplot(2,2,4)
%             plot(29:31,std(ch_Slave(:,29:31)),'-o');
            
        end
        
        pause(0.5);
    end
    
    
elseif strcmp(plot_conf,'sg cal')
    factor_Master=ones(1,length(factor_Master)); %[in volts]
    factor_Slave=ones(1,length(factor_Slave));
    
    for k=1:1000
        ch_multiplex=multiplex_read(path_Master,byts_to_plot);
        ch_Master=multiplex_decoding(ch_multiplex,norm_coeff_Master,ch_offset_Master,factor_Master,jump);
        ch_multiplex=multiplex_read(path_Slave,byts_to_plot);
        ch_Slave=multiplex_decoding(ch_multiplex,norm_coeff_Slave,ch_offset_Slave,factor_Slave,jump);
        
        % ----------S.G Calibration Plot
        % --------acq132_093
        subplot(4,2,1);
        plot(ch_Master(:,[2,3,5,6,8,9,11,12,14,15]),'.-'); %ch with positive  offset
        title('First card-093 - ch 2,3,5,6,8,9,11,12,14,15')
        subplot(4,2,2);
        plot(ch_Master(:,1:3:13),'.-');%ch with negative offset
        subplot(4,2,3);
        plot(ch_Master(:,[18,19,21,22,24,25,27,28,30,31]),'.-');%ch with positive  offset
        title('Second card-093 - ch 18,19,21,22,24,25,27,28,30,31')
        subplot(4,2,4);
        plot(ch_Master(:,17:3:29),'.-');%ch with negative offset
        %
        %---------first card _094
        subplot(4,2,5);
        plot(ch_Slave(:,[2,3,5,6,8,9,11,12,14,15]),'.-'); %ch with positive offset
        title('First card-094 - ch 2,3,5,6,8,9,11,12,14,15')
        subplot(4,2,6);
        plot(ch_Slave(:,1:3:13),'.-');%ch with negative offset
        subplot(4,2,7);
        plot(ch_Slave(:,[18,19,21,22,24,25,27,28]),'.-');%ch with positive  offset
        title('Second card-094 - ch 18,19,21,22,24,25,27,28')
        subplot(4,2,8);
        plot(ch_Slave(:,17:3:28),'.-');%ch with negative offset
        
% plot([1:15 17:31],std(ch_093(:,[1:15 17:31])),'-o');
% hold all
% plot([1:15 17:28],std(ch_094(:,[1:15 17:28])),'-o');
% hold off;
% ylim([0.008 0.014]);
        %    ch1_std=std(ch_094(:,23))
        %    ch2_std=std(ch_094(:,25))
        % %
        pause(0.5);
    end
    elseif strcmp(plot_conf,'elsa sg cal')
    factor_Master=ones(1,length(factor_Master)); %[in volts]
    factor_Slave=ones(1,length(factor_Slave));
    
    for k=1:1000
        ch_multiplex=multiplex_read(path_Master,byts_to_plot);
        ch_Master=multiplex_decoding(ch_multiplex,norm_coeff_Master,ch_offset_Master,factor_Master,jump);
        ch_multiplex=multiplex_read(path_Slave,byts_to_plot);
        ch_Slave=multiplex_decoding(ch_multiplex,norm_coeff_Slave,ch_offset_Slave,factor_Slave,jump);
        
        % ----------S.G Calibration Plot
        % --------acq132_093
        subplot(2,2,1);
        plot(ch_Master(:,1:15),'.-'); %ch with positive  offset
        title('ch Master 1-15')
        subplot(2,2,3);
        plot(ch_Master(:,17:31),'.-');
        title('ch Master 17-31')
        subplot(2,2,2);
        plot(ch_Slave(:,1:15),'.-');%ch with negative offset
        title('ch Slave 1-15')
        subplot(2,2,4);
        plot(ch_Slave(:,17:28),'.-');
        title('ch Slave 17-28')
     
        
% plot([1:15 17:31],std(ch_093(:,[1:15 17:31])),'-o');
% hold all
% plot([1:15 17:28],std(ch_094(:,[1:15 17:28])),'-o');
% hold off;
% ylim([0.008 0.014]);
        %    ch1_std=std(ch_094(:,23))
        %    ch2_std=std(ch_094(:,25))
        % %
        pause(0.5);
    end
end


%--------------format file read
function [norm_coeff]=format_read(path_format)
DELIMITER = '\t';
HEADERLINES = 42;
newData1 = importdata(path_format, DELIMITER, HEADERLINES);
vars = fieldnames(newData1);
for i = 1:length(vars)
    vars{i}= newData1.(vars{i});
end
norm_coeff=vars{1,1};

%----------------read multiplexed data
function [ch_multiplex]=multiplex_read(path,byts_to_plot,origin)
if nargin<3 || strcmp(origin,'eof')
    origin='eof';
    offset=-byts_to_plot;
else
    offset=0;
end

fid= fopen(path,'r');
fseek(fid,offset,origin);
ch_multiplex= fread(fid,byts_to_plot/2,'int16'); %factor 1/2 is needed because 'int16' is 2 bytes
fclose(fid);
%---------------- decode multiplexed data
function [ch]=multiplex_decoding(ch_multiplex,norm_coeff,ch_offset,factor_09,jump)

for j=2:33
    ch(:,j-1)=ch_multiplex(j:33:end);
end
ch1_16=ch(1:jump:end,1:2:end);
ch17_32=ch(1:jump:end,2:2:end);
ch=[ch1_16 ch17_32];
for j=1:32
    ch(:,j)=(ch(:,j)*norm_coeff(j,1)+norm_coeff(j,2)-ch_offset(j))/factor_09(j);
end

% %------Read Ph line
% function line=phantomReadLineEOF()
%
% f = fopen('c:\Frics\tmp\Ph\0.bin','r');
% fseek(f,(-1)*imSize*2,'eof');
% Im=fread(f,imSize,'uint16');
% fclose(f);
% Imr=reshape(Im',PhMeta.ImageWidth,PhMeta.ImageHeight)';
% line=mean(Imr,1);

