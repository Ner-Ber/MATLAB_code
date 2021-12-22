function plot_stress_distribution_at_t(s,t,baseFig)

baseFigStart=baseFig;

spl=calc_Fx_stress_splined(s);
[~,index]=min(abs(s.t-t));

figure(baseFig);plot(spl.x(2:end-1),spl.Fx(index,:),'.-');
baseFig=baseFig+1;
title ('Fx')
figure(baseFig);plot(spl.x,spl.Syy(index,:),'.-');
baseFig=baseFig+1;
title ('Syy')
figure(baseFig);plot(spl.x(2:end-1),spl.Fx(index,:)./spl.Syy(index,2:end-1),'.-');
title ('Fx/Syy')

for j=baseFigStart:baseFig
    hold (get(j,'CurrentAxes'),'all')
end
