function anls_taur_Check(anl)

%x_sgForPlot=[105.36, 113.02, 121.02, 128.4];
%x_sgForPlot=[85 100, 126.28, 149.56];
%x_sgForPlot=[81.5 97 105.5 130];
x_sgForPlot=[89 113 135 ];
%index_e=[ 13, 14, 15];
index_e=1:length(anl.e);
% ------ Compute Uxyf and Uxy0

% for j=1:length(anl.y3_5.x_sg)
%     [~,index]=min(abs(anl.pre.x_sg-anl.y3_5.x_sg(j)));
%
%     anl.y3_5.Uxy0(:,j)=anl.pre.Uxy(:,index);
% end

anl.y3_5.Uxyf=anl.y3_5.Uxy0 - anl.y3_5.Uxy0_f;


%%-----Go over all the events
figN1=figure;
figN2=figure;
for j=1:length(index_e)
    
    e=anls_front_in_space(anl.lgnd{j}(12:19),anl.e(index_e(j)),80,170);
    figure(figN1);
    for k=1:length(x_sgForPlot)
        [~,index_sg1]=min(abs(e.acqE.x_sg-x_sgForPlot(k)));
        [~,index_sg2]=min(abs(anl.y3_5.x_sg-x_sgForPlot(k)));
        subplot(2,2,k);
        plot(e.acqE.intVdtMat(:,index_sg1),e.acqE.Uxy(:,index_sg1)-anl.y3_5.Uxyf(j,index_sg2),'.-');
        %plot(e.acqE.intVdtMat(:,index_sg1),e.acqE.Uxy(:,index_sg1)-anl.y3_5.Uxy0(j,index_sg2),'.-');
        xlim([-50 50]);
        my_legend_add(num2str(e.acqE.event));
        title(num2str(e.acqE.x_sg(index_sg1)));
        hold all;
        pause(0.2);
    end
    
    figure(figN2);
    for k=1:length(x_sgForPlot)
        [~,index_sg1]=min(abs(e.acqE.x_sg-x_sgForPlot(k)));
        [~,index_sg2]=min(abs(anl.y3_5.x_sg-x_sgForPlot(k)));
        subplot(2,2,k);
        tau=calc_dynamic_stress(e.acqE.Uxy(:,index_sg1),anl.y3_5.Uxy0(index_e(j),index_sg2));
        plot(e.acqE.intVdtMat(:,index_sg1),tau,'.-');
        %plot(e.acqE.intVdtMat(:,index_sg1),e.acqE.Uxy(:,index_sg1)-anl.y3_5.Uxy0(j,index_sg2),'.-');
        xlim([-50 50]);
        my_legend_add(num2str(e.acqE.event));
        title(num2str(e.acqE.x_sg(index_sg1)));
        hold all;
        pause(0.2);
    end
    
end

