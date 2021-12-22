function phantomImMeshMovie(exp_dir,eventNum,Tstart,Tinterval,Tend,timeBase,smt,lineNum)
%phantomImMeshMovie(exp_dir,eventNum,Tstart,Tinterval,Tend,timeBase,smt,lineNum)
%eventNum=0 is slow acquisition. [msec] timeBase=1000.
%The function uses phantomGetIms -> reads all ims at onces.


Ims=phantomGetIms(exp_dir,eventNum,Tstart,Tinterval,Tend,timeBase,smt,lineNum);
Ims0=phantomGetIms(exp_dir,eventNum,'start','min',0,timeBase,smt,lineNum);
firstImr=mean(Ims0.ims(:,:,1:20),3);
%Ims.ims=Ims.ims-repmat(firstImr,[1,1,length(Ims.ims(1,1,:))]);

ind=1;

figBase=gcf;
%plotbrowser off;
%close(figBase);


while ind>0&&ind<length(Ims.t)
    
    figure(figBase);
    imagesc(Ims.x,(lineNum-1)/8*5.5,Ims.ims(:,:,ind)./firstImr);
    %imagesc(Ims.x,(lineNum-1)/8*50.5,Ims.ims(:,:,ind));
    set(gca,'ydir','normal')
    axis image
    caxis([0.8 1.02]);
   % my_mesh(Ims.x,lineNum,Ims.ims(:,:,ind),0);
   %caxis([0 400]);
  % caxis([0 2*mean(mean(Ims.ims(:,:,ind)))]);
    
    %axis tight;
    colorbar;
    title(num2str(Ims.t(ind)));
    view(0,90);
    
   %------Comment those lines for for continium movie
   set(gca,'NextPlot','replaceChildren');
   set(0,'DefaultFigureWindowStyle','normal')
   choice = menu('','R','F','FX20','RX20' , 'Quit');
   set(0,'DefaultFigureWindowStyle','docked')
    
%   pause (0.05)
 % choice=2;
    
     if(choice==1)
        ind=ind-1;
    elseif(choice==4)
        ind=ind-20;
    elseif(choice==2)
        ind=ind+1;
    elseif(choice==3)
        ind=ind+20;
    else
        break
    end
    
end



