function anl=SimFindBumpLocation(fault)

% path='C:\Users\IlyaS\Desktop\2015-02-11';
% xMat=csvread([path '\x_fault.csv']);
% x=xMat(1,:);
% t=csvread([path '\time_fault.csv']);
% UxyMat=csvread([path '\eps12_fault.csv']);

x=fault.x(1:end);%x=fault.x(2400:end);
t=fault.t;
UxyMat=fault.Sxy(:,:)';%UxyMat=fault.Sxy(2400:end,:);
%-------Find Bump Location
k=1;
%figure;
%4.0-170,4.05-110-172,4.1-175,4.2-108,4.3-75
%Bilateral-4.2 133:171
%new 4.0-1900 4.05-1200
%NoHealing 4.3-280:1500, 4.2-460-1500, 4.1 840-1700

for j=50:1:length(UxyMat(1,:))
    
    Uxy=UxyMat(:,j);%Uxy should column oriented
    
    [~,indexTip]=max(Uxy);
    xTip(k)=x(indexTip);
    tTip(k)=t(k);
    indexR=find(Uxy>3.704e6,1);%find end of the cohesive zone
    xCend(k)=x(indexR-1)+( x(indexR)-x(indexR-1) )*(1-( Uxy(indexR)-Uxy(indexR-1))/0.5e6);%some interpolation for better resolution
    
%     index=indexTip+40; %Need to play with that number to get good results
%     
%   %  [~,indexBump]=min(abs(Uxy(index:end)-Uxy(end) -50)); %Finde first arrival
%    % indexBump=indexBump+index-1;
%     %xBumpP(k)=x(indexBump);
%     
%     
%    %[indexMax ~]=localMaximum(Uxy(index:end),200); %Need to play with that number to get good results %Uxy should column oriented
%     [c indexMax]=max(Uxy(index:end)); %Need to play with that number to get good results %Uxy should column oriented
%     
%     if( ~isempty(indexMax) )
%         indexBump=indexMax(1)+index-1;
%         
%         xBump(k)=x(indexBump);
%         tBump(k)=t(j);
%         SxyB(k)=Uxy(indexBump);
%         
% %         [~,indexBump]=min(abs(Uxy(index:end)-Uxy(end) -100)); %Finde first arrival
% %         indexBump=indexBump+index-1;
% %         xBumpP(k)=x(indexBump);
%         
%         %         %plot(Uxy,'.-');
%         %     plot(indexBump(j),Uxy(indexBump(j)),'o');
%         %plot(x- xCend(k),Uxy,'.-');
%         
%        %  plot(x- x(indexBump),Uxy);
%         % hold all;
%          
        k=k+1;
  %  end
    
    end

    
%anl.xBump=xBump;
%anl.xBumpP=xBumpP;
%anl.t=tBump;
anl.tTip=tTip;
%anl.SxyB=SxyB;
anl.xTip=xTip;
anl.xCend=xCend;
anl.v=diff(xTip)./diff(anl.tTip);


%Find Anti_Bump
% [~,indexX]=min(abs(fault.x-0.0585));
% [~,indexT]=min(abs(fault.t-1.02e-4));
% [~,indexT2]=min(abs(fault.t-2e-4));
%for j=1:(indexT2-indexT)  [S(j),index(j)]=min(fault.Sxy(indexX:end,j));index(j)=index(j)+indexX-1;end


% %-------read data
% path='C:\Users\IlyaS\Desktop\2015-02-11\mesh';
%
% x=csvread([path '\x_map.csv'],0,0,[0 0 6300 0]);
% y=csvread([path '\y_map.csv'],0,3200,[0 3200 0 3200])
% Uxy=csvread([path '\eps12_map.csv'],0,0,[0 0 6300 3200]);
%






