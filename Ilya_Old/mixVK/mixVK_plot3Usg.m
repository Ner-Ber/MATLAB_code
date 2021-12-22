function mixVK_plot3Usg(exp,event,esg,V,K)

event_vec=event;
current_dir=pwd;
lastSlashPosition = find(current_dir == '\', 1, 'last');
Date = current_dir(lastSlashPosition+1:end);

figure;hold all;
for k=1:length(event_vec)
    event=event_vec(k);
    e=esg;
%     e=e_exp.e{event};
    plot(e.acqEc.intVdtMat(:,V),e.acqEc.Uxx(:,V)+e.acqEc.Uxy(:,V)-repmat(e.pre.Uxx(V)+e.pre.Uxy(V),length(e.acqEc.Uyy),1),'x-','MarkerSize',3)
    plot(e.acqEc.intVdtMat(:,K),e.acqEc.Uxx(:,K)+e.acqEc.Uxy(:,K)-repmat(e.pre.Uxx(K)+e.pre.Uxy(K),length(e.acqEc.Uyy),1),'.-','MarkerSize',6)
    legend([num2str([e.acqE.x_sg(V) e.acqE.x_sg(K)]')]); legend off
    title([Date ' ' exp ' ' num2str(event) ' Uyy'],'FontSize',14,'FontWeight','bold')
    xlim([-60 60])
end

figure;hold all;
for k=1:length(event_vec)
    event=event_vec(k);
    e=esg;
%     e=e_exp.e{event};
    plot(e.acqEc.intVdtMat(:,V),e.acqEc.Uxx(:,V)-repmat(e.pre.Uxx(V),length(e.acqEc.Uyy),1),'x-','MarkerSize',3)
    plot(e.acqEc.intVdtMat(:,K),e.acqEc.Uxy(:,K)-repmat(e.pre.Uxy(K),length(e.acqEc.Uyy),1),'.-','MarkerSize',6)
    legend([num2str([e.acqE.x_sg(V) e.acqE.x_sg(K)]')]); legend off
    title([Date ' ' exp ' ' num2str(event) ' Uxx'],'FontSize',14,'FontWeight','bold')
    xlim([-60 60])

end

figure;hold all;
for k=1:length(event_vec)
    event=event_vec(k);
    e=esg;
%     e=e_exp.e{event};
    plot(e.acqEc.intVdtMat(:,V),e.acqEc.Uxy(:,V)-repmat(e.pre.Uxy(V),length(e.acqEc.Uyy),1),'x-','MarkerSize',3)
    plot(e.acqEc.intVdtMat(:,K),e.acqEc.Uxx(:,K)-repmat(e.pre.Uxx(K),length(e.acqEc.Uyy),1),'.-','MarkerSize',6)
    legend([num2str([e.acqE.x_sg(V) e.acqE.x_sg(K)]')]); legend off
    title([Date ' ' exp ' ' num2str(event) ' Uxy'],'FontSize',14,'FontWeight','bold')
    xlim([-60 60])

end