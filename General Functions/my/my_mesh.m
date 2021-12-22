function my_mesh(x,t,D,normalize)
%normalize=1 -> FirstLine is subtructed:               D=D-FirstLine;
%normalize=2 -> The data normolized by FirstLine:      D=D./FirstLine;
%else The original data is plotted

if(length(t(1,:))>1 && length(t(:,1))==1)
    t=t';
end
% 
% x=[x x(end)+(x(end)-x(end-1))]; %mesh makes troubles with the last column on plot
% D=[D D(:,end)];
% 
% if (length(t)>2)
% t=[t ; t(end)+(t(end)-t(end-1))]; %mesh makes troubles with the last column on plot
% D=[D ; D(end,:)];
% end

if (nargin<4)
    normalize=0;
end

if normalize==1 %FirstLine is subtructed
    FirstLine=repmat(mean(D(1:50,:)),length(t),1);
    D=D-FirstLine;
else if normalize==2%The data normolized by FirstLine
        FirstLine=repmat(mean(D(1:50,:)),length(t),1);
        D=D./FirstLine;
    end
end


%------------------meshes
%  X=repmat(x,length(t),1);
%  T=repmat(t,1,length(x));
% % %caxis
%  mesh(X,T,D,...
%     'FaceColor','interp',...
%     'EdgeColor','interp');view(0,90);
%

%mesh(X,T,D,...
   % 'FaceColor','flat');view(0,90);
%uimagesc(x,t,D);
imagesc(x,t,D);
%imagesc(x,t,(D*1280-320)/(1280-320));
set(gca,'ydir','normal')
colormap('jet')

axis tight;
colorbar
