function mu=calc_mu(acqE)

Syy=mean(acqE.Syy(1:500,:),1);
mu=acqE.Sxy./repmat(Syy,length(acqE.t),1);