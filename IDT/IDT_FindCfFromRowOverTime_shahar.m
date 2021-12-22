function [f1,f2,FrontTime,StairLoc,FrontVel,VelLoc]=IDT_FindCfFromRowOverTime_shahar(RowOverTime,TimeCrop,SpatialCrop,smt)

framerate=200000;
pixelsize=86.5*1e-6;


if (nargin<4)                                    %% default is smooth factor 4
    smt=4;
end

if (nargin<3 || strcmp(SpatialCrop,'all'))       %% default is check for all the spatial pixels
    SpatialCrop=1:size(RowOverTime,2);
end

if (nargin<2 || strcmp(TimeCrop,'all'))          %% default is check for all frames
    TimeCrop=1:size(RowOverTime,1);
elseif strcmp(TimeCrop,'EventArea')              %% define the TimeCrop as the close time around the event
    M=mean(RowOverTime,2);
    EventFrame=find(M==min(M));
    frames=100;                                  %% take an area of 100 frames
    TimeCrop=EventFrame-frames:EventFrame+frames;
end

%--smooth the picture
SmoothedROT=conv2(RowOverTime,ones(smt)/smt^2,'same');

%--Crop Image
CroppedROT=SmoothedROT(TimeCrop,SpatialCrop);

%--2nd derivative
ScndDrvROT=conv2(CroppedROT,[1 -2 1]','valid');


f1=[];
f2=[];
% FrontLoc=[];
% problematic=[];
% k=1;

% n = 4; %remove from edge
% [~,Imxs]=max(ScndDrvROT(n:end-n,:),[],1);
% [~,Imns]=min(ScndDrvROT(n:end-n,:),[],1);
% Imxs = Imxs+n-1;
% Imns = Imns+n-1;

% colmn=ScndDrvROT(:,I);
% f1(i)=find(colmn==maxs,1);
% f2(i)=find(colmn==mins,1,'last');

for i=1:size(CroppedROT,2)
    maxs=max(ScndDrvROT(:,i));
    mins=min(ScndDrvROT(:,i));
    colmn=ScndDrvROT(:,i);
    f1(i)=find(colmn==maxs,1);
    f2(i)=find(colmn==mins,1,'last');
    
    %     yy=row(f2(i):f1(i));
    %     xx=f2(i):f1(i);
    %     try
    %         FrontLoc(i)=interp1(yy,xx,0);
    %     catch
    %         FrontLoc(i)=FrontLoc(i-1);
    %         problematic(k)=i;
    %         k=k+1;
    %     end
    
    
end

%--find midpoint
fronLoc=floor((Imxs+Imns)./2);    %this produces a "stairs" graph

%--take only the startpoint of the stair
Stairs=find(diff(fronLoc));
FrontTime=fronLoc(Stairs);
% SpatVec=1:size(RowOverTime,2);
StairLoc=Stairs+SpatialCrop(1);


%--Units
StairLoc=StairLoc.*pixelsize;
FrontTime=FrontTime./framerate;

%--compute velocity
dX=diff(StairLoc);
dT=diff(FrontTime);
FrontVel=dX./dT;
VelLoc=movmean(StairLoc,2,'Endpoints','discard');

figure
subplot(2,1,1);imagesc(RowOverTime)
subplot(2,1,2);plot(VelLoc*1e3,FrontVel,'x')
xlim([1 size(RowOverTime,2)].*pixelsize*1e3)
ylabel('Velocity [m/s]')
xlabel('Location [mm]')


end