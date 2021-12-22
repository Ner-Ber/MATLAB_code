function outmu=Odsv_2mu(path,involts)
%First it loads calibration parameters V, d, p for example from ex_dir\cal\file, interpolates
%it and returns the corrcted outmu

%-------Method to create calibration mfile. calibrate acqS.ch_094_ex(:,4) with acqS.ch_094_ex(:,2)
%0.use index=logical((acqS.ch_094_ex(:,4)<2) .*(acqS.ch_094_ex(:,4)>0.21)); to define upper and lower limits.
%1. plot(acqS.ch_094_ex(index,4),acqS.ch_094_ex(index,2));
%2. fig=myget_axis();
%3. define xx 0.2:0.01:2
%4.use csaps to interpolate->  yy=csaps(fig.x,fig.y,p,xx);
%5. hold all; plot(xx,yy); until fits well change p.
%6. V=fig.x; d=fig.y;
%7.save the output parameters (V,d,p) at exp_dir\cal\Odscal_ch_094_ex4 (for external ch4) 

%----------------
%cal=load([exp_dir '\cal\' filename]);
cal=load(path);
outmu=csaps(cal.V,cal.d,cal.p,involts);

