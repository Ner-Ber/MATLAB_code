function [D]=subtruct_norm(D,normalize,v)
%The function subtructs inital values from the data or normolizes it
%according initial values.
%If v entered, it is used insread of the initial values.  
%normalize=1 -> The data normolized by FirstLine:      D=D./FirstLine
%else subtructs FirstLine.

if nargin<2
    normalize=0;
end

if nargin~=3
    index_end=5;%number of points to get the mean value for substruction and normalization
    v=mean(D(1:index_end,:),1);
elseif length(v(1,:))==1 
    v=v';
end


%[row,col]=size(D);
% if row==1 || col==1
%     
%     if row==1
%         %one dimension data
%         FirstLine=repmat(mean(D(1:index_end)),1,length(D));
%     else
%         FirstLine=repmat(mean(D(1:index_end)),length(D),1);
%     end
%     
%     if normalize == 1
%         D=D./FirstLine;
%     else
%         D=D-FirstLine;
%     end
%     
% else
    % Matrice data form
    
    FirstLine=repmat(v,length(D(:,1)),1);  %mean treats the column as vectors
    if normalize == 1
        D=D./FirstLine;
    else
        D=D-FirstLine;
    end
