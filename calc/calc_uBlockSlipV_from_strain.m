function spl=calc_uBlockSlipV_from_strain(struct,frontV,UxySmt)
%frontV is row vector. It's length should equal the number of strain gages

if (nargin<3)
UxySmt=floor(71./frontV*157);
end

UxxSmt=1;

%------------Read exp_details.txt
expDetails=expDetailsRead(struct.exp);
height=expDetails.sg_height;%[mm]
%---------------

Uxy=my_smooth(struct.Uxy,UxySmt);
Uxx=my_smooth(subtruct_norm(struct.Uxx),UxxSmt);
t=struct.t;
clear struct;

fudgeFactor=1;
spl.uBlockSlipV1=-repmat(frontV,length(t)-1,1).*Uxx(1:end-1,:)*100/1000; %[V]=m/sec -> cm/sec, [Uxx]=mstrain -> strain
spl.uBlockSlipV2=+repmat(height,length(t)-1,1)*fudgeFactor.*diff(Uxy,1,1)./repmat(diff(t),1,length(Uxy(1,:)))*100/1000/1000*1000; %100-> cm/sec ,[height]=mm ,[Uxy]= mstrain [t]=msec, %two '-' 
spl.uBlockSlipV=spl.uBlockSlipV1+spl.uBlockSlipV2;
spl.t=t(1:end-1);


%plot(spl.t,[spl.uBlockSlipV1(:,6) spl.uBlockSlipV2(:,6) spl.uBlockSlipV(:,6)],'.-');