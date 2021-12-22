function e=anls_single_location(exp,event,Tstart,Tend,sg_num,frontV,baseFig)
%if baseFig not entered -> no plot
smtUblockV=1;
smtuBlockVFromStrain=1;
smtFx=3;

frontV=ones(1,15)*frontV;

e=acq132_event_get_data(exp,event,Tstart,Tend,1,'Sxy','Syy','Sxx','x_sg','ch_094_ex','Uxx','Uxy','Uyy');
% e.dBlock=smooth(e.ch_094_ex(:,3),1);
% e.uBlock=smooth(e.ch_094_ex(:,5),1);

% e.dBlock=smooth(e.ch_094_ex(:,2),1);
% e.uBlock=smooth(e.ch_094_ex(:,4),1);

e.dBlock=smooth(e.ch_094_ex(:,2),1);
e.uBlock=smooth(e.ch_094_ex(:,5),1);

slip=(e.uBlock-e.dBlock);
velocity=diff(smooth(slip,1))./diff(e.t)/10;% [cm /sec]
e.slip=slip;
e.v=velocity;
e.t=(e.t);
e.dBlockV=diff(smooth(e.dBlock,1))./diff(e.t)/10;% [cm /sec]
e.uBlockV=diff(smooth(e.uBlock,smtUblockV))./diff(e.t)/10;% [cm /sec]

e.uBlockVFromStrain=calc_uBlockSlipV_from_strain(e,frontV,smtuBlockVFromStrain);
e.uBlockSlipFromStrain=cumsum(e.uBlockVFromStrain.uBlockSlipV(:,sg_num))/10^6*10^4;
e.uBlockVFromStrainNeg=calc_uBlockSlipV_from_strain(e,-frontV,1);

e.Fx=calc_Fx_by_time_deriv(e,frontV,smtFx);
e.FxNeg=calc_Fx_by_time_deriv(e,-frontV,smtFx);

if  exist([exp '\Ph'],'dir')==7
    phE=phantomGetLines(exp,event,'start','min','end',1000,1);
    e.A=get_A_at_x(phE,e.x_sg,11);
end

% if  exist([exp '\vds'],'dir')==7
%     vds=vds_event_get_data(exp,event,100,400,3);
%     figure(baseFig);
%     my_mesh(vds.x,vds.t,vds.lines,2);
%     baseFig=baseFig+1;
% end

if nargin==7;
    figure(baseFig);
    baseFig=baseFig+1;
    plotyy(e.t,e.uBlock,e.t,subtruct_norm([e.Uxx(:,sg_num) e.Uxy(:,sg_num) e.Uyy(:,sg_num)]))  
   
   figure(baseFig);
   baseFig=baseFig+1;
   plotyy(e.A.t,subtruct_norm(e.A.lines(:,sg_num),1),e.Fx.t,subtruct_norm([e.Fx.Fx(:,sg_num) e.Fx.Sxx(:,sg_num) e.Fx.Sxy(:,sg_num) e.Fx.dSxxdx(:,sg_num)])) 
   title('A Fx Sxx Sxy dSxxdx');
   
   figure(baseFig);
    baseFig=baseFig+1;
    plot(e.t(2:end),e.uBlockV,'.-',e.uBlockVFromStrain.t,[e.uBlockVFromStrain.uBlockSlipV(:,sg_num) e.uBlockVFromStrainNeg.uBlockSlipV(:,sg_num)],'.-')  
   
   figure(baseFig);
   baseFig=baseFig+1;
   ind=7:14;
   plot(e.Fx.t,subtruct_norm([e.Fx.Fx(:,ind) e.Fx.Sxx(:,ind) e.Fx.Sxy(:,ind) e.Fx.dSxxdx(:,ind)]),'.-') 
   legend({num2str([e.x_sg(ind) e.x_sg(ind) e.x_sg(ind) e.x_sg(ind)]')});legend('off');
   title('Fx Sxx Sxy dSxxdx');
   
   figure(baseFig);
   baseFig=baseFig+1;
   plotyy(e.A.t,subtruct_norm(e.A.lines(:,:),1),e.t,subtruct_norm([e.Uxy(:,:) e.Uxx(:,:)])) 
  % plotyy(e.A.t,subtruct_norm(e.A.lines(:,[7 9:end]),1),e.Fx.t,subtruct_norm([e.Fx.Fx(:,[7 9:end]) e.Fx.Sxx(:,[7 9:end])])) 
   legend({num2str([e.x_sg e.x_sg e.x_sg]')});legend('off');
   title('A Uxy Uxx');
end




