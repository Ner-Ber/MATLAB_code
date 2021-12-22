function dtip=data_tip_join(dtip1,dtip2)

numrow1=size(dtip1.x,1);
numcol1=size(dtip1.x,2);
numrow2=size(dtip2.x,1);
numcol2=size(dtip2.x,2);
%create full new structure
dtip.x=[dtip1.x ;NaN(numrow2,numcol1)];
dtip.y=[dtip1.y ;NaN(numrow2,numcol1)];

%join intersected events
[~, idtip1, idtip2] = intersect(dtip1.event,dtip2.event);
dtip.x(numrow1+1:end,idtip1)=dtip2.x(:,idtip2);
dtip.y(numrow1+1:end,idtip1)=dtip2.y(:,idtip2);

%sort intersected events with respect to x
[dtip.x(:,idtip1),s]=sort(dtip.x(:,idtip1));

for j=1:length(s(1,:))
dtip.y(:,idtip1(j))=dtip.y(s(:,j),idtip1(j));
end

%clear rows full of NaN
for j=numrow1+numrow2:-1:1
    ind=isnan(dtip.x(j,:));
    if nnz(ind)~=numcol1; %columns are sorted,i.e, only at the end of the matrix posibble NaN row
        break
    else
        dtip.x(j,:)=[];
        dtip.y(j,:)=[];
    end
    
end

%------join the rest of the events 
%find not intersected events
idtip22=idtip2;
idtip2=1:length(dtip2.x(1,:));
idtip2(idtip22)=[];%index vec of different events
%join  NaN matrix
dtip.x=[dtip.x NaN(size(dtip.x,1),length(idtip2))];
dtip.y=[dtip.y NaN(size(dtip.x,1),length(idtip2))];
%add the non intersected events at the end of the matrix
dtip.x(:,numcol1+1:end)=dtip2.x(:,idtip2);
dtip.y(:,numcol1+1:end)=dtip2.y(:,idtip2);
dtip.event=[dtip1.event dtip2.event(idtip2)];
%sort according to event number
[dtip.event,s]=sort(dtip.event);
dtip.x=dtip.x(:,s);
dtip.y=dtip.y(:,s);







