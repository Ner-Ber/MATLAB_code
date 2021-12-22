function out=calc_uBlockSlipV_from_strain_Integrated(struct,tStart,tEnd,smt)


[~,indexStart]=min(abs(struct.t-tStart));
[~,indexEnd]=min(abs(struct.t-tEnd));
out.t=struct.t(indexStart:indexEnd);

struct.Uxx=my_smooth(subtruct_norm(struct.Uxx),smt);
struct.Uxx=-(struct.Uxx-repmat(mean(struct.Uxx(1000:2000,:),1),length(struct.Uxx(:,1)),1));
Uxx=struct.Uxx(:,1:end);
x_sg=struct.x_sg(1:end-2);

%find shift and correction betwean sequent sg.
for j=1:length(x_sg)-1
    [c,lags]=xcorr(Uxx(indexStart:indexEnd,j),Uxx(indexStart:indexEnd,j+1),'coeff');%It's better to cut before xcorr
    [~,index]=max(c);
    shift(j)=lags(index);
    Uxx_shifted=(circshift(struct.Uxx(:,j+1),shift(j)));
    Uxx_shifted1=(circshift(struct.Uxx(:,j+1),sum(shift(1:j))));
    deltaUxx(:,j)=(Uxx_shifted-Uxx(:,j))/abs(shift(j));
    figure(10);plot(struct.t,[Uxx(:,j), Uxx_shifted],'.-');hold all;
    %figure(10);plot(struct.t,Uxx_shifted1,'.-');hold all;
end

%deltaUxx=deltaUxx*0;
shift=-shift;
UxxInterp=ones(length(Uxx(:,1)),sum(abs(shift))+1);
x=ones(1,sum(abs(shift))+1);
deltaX=diff(x_sg)./abs(shift);

%Create interpolation for Uxx betwean all two sg.
index=0;
for j=1:length(x_sg)-1
    index=index+1;
    x(index)=x_sg(j);
    UxxInterp(:,index)=Uxx(:,j);
    UxxTmp=UxxInterp(:,index);
    
    for k=1:abs(shift(j))-1
        index=index+1;
        UxxTmp=UxxTmp+deltaUxx(:,j);
        UxxInterp(:,index)=UxxTmp;
        UxxInterp(:,index)=circshift(UxxInterp(:,index),k);
        x(index)=x(index-1)+deltaX(j);
    end
end
UxxInterp(:,index+1)=Uxx(:,j+1);
x(index+1)=x_sg(j+1);
out.x=x;

%create slip and slip velocity
for j=1:length(x_sg)-1
    [~,index]=min(abs(x-x_sg(j)));
    out.Ux(:,j)=trapz(x(index:end),UxxInterp(indexStart:indexEnd,index:end),2);%[x]=mm,[Uxx]=mstrain->[Ux]=mu
    
end
out.Vx=diff(out.Ux,1)*100'; %[Vx]= cm/sec
out.UxxInterp=UxxInterp(indexStart:indexEnd,:);



