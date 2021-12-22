function anls_Af_Check(anl)

x_sgForPlot=[105.36, 121.02, 130, 136];
%x_sgForPlot=[85 100, 126.28, 149.56];
%x_sgForPlot=[81.5 97 105.5 130];
%x_sgForPlot=[89 113 135 ];
%index_e=[ 13, 14, 15];
anl.e=anl.event;
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

for j=1:length(index_e)
    
    e=anls_front_in_space(anl.lgnd{j}(12:19),anl.e(index_e(j)),80,170);
    figure(figN1);
    
    for k=1:length(x_sgForPlot)
        
        [~,index_sg2]=min(abs(anl.y3_5.x_sg-x_sgForPlot(k)));
        
        xStart=anl.y3_5.x_sg(index_sg2)-5;
        xEnd=anl.y3_5.x_sg(index_sg2)+5;
        [~,index1]=min(abs( e.phECut.frontX-xStart ));
        %[~,index2]=min(abs( e.phECut.frontX-xEnd ));
        index2=index1;
        
        subplot(2,2,k);
        
        for l=index1:index2
            plot(e.phECut.xOffset(l,:),(e.phECut.lines(l,:)-e.phECut.Af(l))/(1-e.phECut.Af(l)),'.-');
            hold all;
            pause(0.2);
        end
        
        xlim([-12 10]);
        ylim([-0.1 1.1]);
        my_legend_add(num2str(e.acqE.event));
        title(num2str(anl.y3_5.x_sg(index_sg2)));
        
        
    end
    
    
end

