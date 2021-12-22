function sum=crosscorr_lines(shift,scan)
    shift=round(shift);
    if(shift>0)
      pos1=scan(1).pos(shift:end);
      depth1=scan(1).depth(1:end-shift+1);
      pos2=scan(2).pos(shift:end);
      depth2=scan(2).depth(1:end-shift+1);
    elseif(shift<0)
      pos1=scan(1).pos(1:end+shift+1);
      depth1=scan(1).depth(-shift:end);
      pos2=scan(2).pos(1:end+shift+1);
      depth2=scan(2).depth(-shift:end);
    else
      pos1=scan(1).pos;
      depth1=scan(1).depth;
      pos2=scan(2).pos;
      depth2=scan(2).depth;      
    end


%zeros are nans
depth1(depth1 == 0)=nan;
depth2(depth2 == 0)=nan;

%clean nans together from both
nans1=isnan(depth1);
depth2(nans1)=[];
depth1(nans1)=[];
pos1(nans1)=[];
pos2(nans1)=[];

nans2=isnan(depth2);
depth2(nans2)=[];
depth1(nans2)=[];
pos1(nans2)=[];
pos2(nans2)=[];

%get unique xs
[x1 inds1]=unique(pos1);
[x2 inds2]=unique(pos2);


%setup vectors for cross corr
d1=depth1(inds1);
d2=interp1(x2,depth2(inds2),x1);

%clean nans together from both
d2(isnan(d1))=[];
d1(isnan(d1))=[];
%
d1(isnan(d2))=[];
d2(isnan(d2))=[];

%now sum square differences of two depths
d=d1-d2;
%sum squares
sum=d'*d/length(d);