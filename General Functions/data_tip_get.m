function f=data_tip_get 
%1.Put cursors on the fig
%3. The function sorts the datatips according to:
%             a. plot curve order b. x axis in each plot curve
%4. OutPut: stractur with x and y matrice. each column -> same plot curve.
%each row -> diferent datatip 

%------usefull commands
%delete(findall(gca,'Type','hggroup','HandleVisibility','off'));
%createDatatip(datacursormode,ax)
%for j=1:length(ax) for k=1:2 createDatatip(datacursormode,ax(j)); end ;end


%---get the figures
line_handles=findobj(gcf,'Type','line');%all the graph handles
line_handles=line_handles(end:-1:1);

%---get the cursors 
dcm_obj = datacursormode(gcf);
cursors = getCursorInfo(dcm_obj);

%---find cursor's line, position and time 
[~,cursor_line]=ismember([cursors.Target],line_handles); 
cP=[cursors.Position];%cP for cursorPosition
cP=reshape(cP,2,length(cP)/2);

%------sort the cursor data according to line (plot curve) order
[cursor_line,s]=sort(cursor_line);
% cP=cP(:,s);

%rearrange the data to matrix,i.e, each column (in
%x and y) is the same line (plot curve), different row is different cursor

index=find(diff(cursor_line)>=1);
index=[1 index+1 length(cursor_line)+1] ;
% numCursorsPerLine=diff(index);
x=NaN(max(numCursorsPerLine),length(index)-1);
y=NaN(max(numCursorsPerLine),length(index)-1);

for j=1:length(index)-1
x(1:numCursorsPerLine(j),j)=cP(1,index(j):index(j)+numCursorsPerLine(j)-1)';
y(1:numCursorsPerLine(j),j)=cP(2,index(j):index(j)+numCursorsPerLine(j)-1)';
end

line_list=cursor_line(index(1:end-1));

%--------sort the data at each plot curve column  according to time
[x,s]=sort(x,1);
for j=1:length(s(1,:))
y(:,j)=y(s(:,j),j);
end

f=struct('line_list',line_list,'x',x,'y',y);
%f=struct('x',x,'y',y);

