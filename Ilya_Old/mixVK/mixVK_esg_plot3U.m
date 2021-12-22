function mixVK_esg_plot3U(exp,event,esg,V,K)

current_dir=pwd;
lastSlashPosition = find(current_dir == '\', 1, 'last');
Date = current_dir(lastSlashPosition+1:end);

% figure;hold all;
% plot(esg.acqE.tOffset(:,V),subtruct_norm(esg.acqE.Uxx(:,V)),'x-','MarkerSize',3)
% plot(esg.acqE.tOffset(:,K),subtruct_norm(esg.acqE.Uxy(:,K)),'.-','MarkerSize',6)
% legend([num2str([esg.acqE.x_sg(V) esg.acqE.x_sg(K)]')]); legend off
% title([Date ' ' exp ' ' num2str(event) ' 1K 3V'],'FontSize',14,'FontWeight','bold')
% xlim([-0.5 0.5])
% 
% 
% figure;hold all;
% plot(esg.acqE.tOffset(:,V),subtruct_norm(esg.acqE.Uxy(:,V)),'x-','MarkerSize',3)
% plot(esg.acqE.tOffset(:,K),subtruct_norm(esg.acqE.Uxx(:,K)),'.-','MarkerSize',6)
% legend([num2str([esg.acqE.x_sg(V) esg.acqE.x_sg(K)]')]); legend off
% title([Date ' ' exp ' ' num2str(event) ' 3K 1V'],'FontSize',14,'FontWeight','bold')
% xlim([-0.5 0.5])
% 
% 
% figure;hold all;
% plot(esg.acqE.tOffset(:,V),subtruct_norm(esg.acqE.Uxy(:,V)+esg.acqE.Uxx(:,V)),'x-','MarkerSize',3)
% plot(esg.acqE.tOffset(:,K),subtruct_norm(esg.acqE.Uxx(:,K)+esg.acqE.Uxy(:,K)),'MarkerSize',6)
% legend([num2str([esg.acqE.x_sg(V) esg.acqE.x_sg(K)]')]); legend off
% title([Date ' ' exp ' ' num2str(event) ' 1+3'],'FontSize',14,'FontWeight','bold')
% xlim([-0.5 0.5])
figure;hold all;
plot(esg.acqE.intVdtMat(:,V),subtruct_norm(esg.acqE.Uxx(:,V)),'x-','MarkerSize',3)
plot(esg.acqE.intVdtMat(:,K),subtruct_norm(esg.acqE.Uxy(:,K)),'.-','MarkerSize',6)
legend([num2str([esg.acqE.x_sg(V) esg.acqE.x_sg(K)]')]); legend off
title([Date ' ' exp ' ' num2str(event) ' 1K 3V'],'FontSize',14,'FontWeight','bold')



figure;hold all;
plot(esg.acqE.intVdtMat(:,V),subtruct_norm(esg.acqE.Uxy(:,V)),'x-','MarkerSize',3)
plot(esg.acqE.intVdtMat(:,K),subtruct_norm(esg.acqE.Uxx(:,K)),'.-','MarkerSize',6)
legend([num2str([esg.acqE.x_sg(V) esg.acqE.x_sg(K)]')]); legend off
title([Date ' ' exp ' ' num2str(event) ' 3K 1V'],'FontSize',14,'FontWeight','bold')



figure;hold all;
plot(esg.acqE.intVdtMat(:,V),subtruct_norm(esg.acqE.Uxy(:,V)+esg.acqE.Uxx(:,V)),'x-','MarkerSize',3)
plot(esg.acqE.intVdtMat(:,K),subtruct_norm(esg.acqE.Uxx(:,K)+esg.acqE.Uxy(:,K)),'MarkerSize',6)
legend([num2str([esg.acqE.x_sg(V) esg.acqE.x_sg(K)]')]); legend off
title([Date ' ' exp ' ' num2str(event) ' 1+3'],'FontSize',14,'FontWeight','bold')

