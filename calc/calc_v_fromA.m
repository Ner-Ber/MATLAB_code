function [frontX frontV ]=calc_v_fromA(exp_dir,eventNum,Tstart,Tend,Xstart,Xend,smt)

phE=phantomGetLines(exp_dir,eventNum,Tstart,'min',Tend,1000,smt,2:4);
x=phE.x;
t=phE.t;
lines=phE.lines;
firstLine=repmat(lines(1,:),length(lines(:,1)),1);

[~,pixStart]=min(abs(phE.x-Xstart));
[~,pixEnd]=min(abs(phE.x-Xend));

A=lines./firstLine;
k=0;

for j=20:length(A(:,1))
    
    %-----find maxima
    % A_tresh_min=1.5e-3;
    % A_tresh_max=1;
    %
    % A_smooth=smooth(A(j,:),smt);
    % diff_A=diff(A_smooth);
    % A_locs=localMaximum(diff_A,50);
    % A_locs_tresh=logical((A_locs>pixStart).*(A_locs<pixEnd).*(diff_A(A_locs)>A_tresh_min).*(diff_A(A_locs)<A_tresh_max));
    % A_locs=A_locs(A_locs_tresh);
    %--------
        
    A_locs=find(A(j,:)>0.95 & A(j,:)<0.97);
    A_locs=A_locs(A_locs>pixStart & A_locs<pixEnd);
    
    if ~isempty(A_locs)
        A_locs=A_locs(ceil(length(A_locs)/2));
        k=k+1;
        xx(k)=x(A_locs);
        tt(k)=t(j);
        
             
% %        
         plot(xx(k),A(j,A_locs),'o',x,A(j,:),'.-');
         %plot(x-xx(k),A(j,:),'.-');
         hold all;
        
    end
end

frontV=diff(xx)./diff(tt);
frontX=xx(1:end-1)+diff(xx)/2;