function mixVK_plot3U(exp,event,esg,choiceV,choiceK,smoothV,smoothK)

if nargin<6
    smoothV=1;
    smoothK=1;
elseif nargin<7
    smoothK=1;
end
event_vec=event;
current_dir=pwd;
lastSlashPosition = find(current_dir == '\', 1, 'last');
Date = current_dir(lastSlashPosition+1:end);

figure;hold all;
for k=1:length(event_vec)
    event=event_vec(k);
    e=esg;
%     e=e_exp.e{event};
    plot(e.acqEc.intVdtMat(:,choiceV),my_smooth(e.acqEc.Uyy(:,choiceV),smoothV)-repmat(e.pre.Uyy(choiceV),length(e.acqEc.Uyy),1),'blackx-','MarkerSize',3)
    plot(e.acqEc.intVdtMat(:,choiceK),my_smooth(e.acqEc.Uyy(:,choiceK),smoothK)-repmat(e.pre.Uyy(choiceK),length(e.acqEc.Uyy),1),'red.-','MarkerSize',6)
    legend([num2str([e.acqE.x_sg(choiceV) e.acqE.x_sg(choiceK)]')]); legend off
    title([Date ' ' exp ' ' num2str(event) ' Uyy'],'FontSize',14,'FontWeight','bold')
    xlim([-60 60])
end

figure;hold all;
for k=1:length(event_vec)
    event=event_vec(k);
    e=esg;
%     e=e_exp.e{event};
    plot(e.acqEc.intVdtMat(:,choiceV),my_smooth(e.acqEc.Uxx(:,choiceV),smoothV)-repmat(e.pre.Uxx(choiceV),length(e.acqEc.Uyy),1),'blackx-','MarkerSize',3)
    plot(e.acqEc.intVdtMat(:,choiceK),my_smooth(e.acqEc.Uxx(:,choiceK),smoothK)-repmat(e.pre.Uxx(choiceK),length(e.acqEc.Uyy),1),'r.-','MarkerSize',6)
    legend([num2str([e.acqE.x_sg(choiceV) e.acqE.x_sg(choiceK)]')]); legend off
    title([Date ' ' exp ' ' num2str(event) ' Uxx'],'FontSize',14,'FontWeight','bold')
    xlim([-60 60])

end

figure;hold all;
for k=1:length(event_vec)
    event=event_vec(k);
    e=esg;
%     e=e_exp.e{event};
    plot(e.acqEc.intVdtMat(:,choiceV),my_smooth(e.acqEc.Uxy(:,choiceV),smoothV)-repmat(e.anl.Uxyf(choiceV),length(e.acqEc.Uyy),1),'blackx-','MarkerSize',3)
    plot(e.acqEc.intVdtMat(:,choiceK),my_smooth(e.acqEc.Uxy(:,choiceK),smoothK)-repmat(e.anl.Uxyf(choiceK),length(e.acqEc.Uyy),1),'r.-','MarkerSize',6)
    legend([num2str([e.acqE.x_sg(choiceV) e.acqE.x_sg(choiceK)]')]); legend off
    title([Date ' ' exp ' ' num2str(event) ' Uxy'],'FontSize',14,'FontWeight','bold')
    xlim([-60 60])

end