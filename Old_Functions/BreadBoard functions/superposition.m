choose_videos=[5 6 7];
MassCenter_cell=cell(4,7);

%% pick the different events
% why do you repeat this 7 times?
% this is for superposition of events from different videos
% for this code you need to change the parameter 'video number' in the...
... 'identify 2D movement' code to be 'i'
A=[];
for i=choose_videos       
    identify_2D_movement;
    MassCenter_cell{1,i}=cell(1,length(coordinatesCell));
    MassCenter_cell{2,i}=cell(1,length(coordinatesCell));
    MassCenter_cell{3,i}=cell(1,length(coordinatesCell));
    MassCenter_cell{4,i}=cell(1,length(coordinatesCell));
    a=[];
    for regionIdx=1:length(coordinatesCell)
        MassCenter_cell{1,i}{regionIdx}=Xcenter(regionIdx,:);
        MassCenter_cell{2,i}{regionIdx}=Ycenter(regionIdx,:);
        MassCenter_cell{3,i}{regionIdx}=diff(Xcenter(regionIdx,:));
        MassCenter_cell{4,i}{regionIdx}=diff(Ycenter(regionIdx,:));
        a(regionIdx)=max(MassCenter_cell{3,i}{regionIdx});
    end
    A=[A a];
end




%% plots
%---the plots are normalized by the mean position of the mass-center...
...before the event. 
%---the plots are synchronized. the event starts at frame 31  
    

%---plot1 - superposition of the movement in X direction
figure
hold all
asparity_number=[];
for i=choose_videos
    for j=1:length(MassCenter_cell{1,i})
        mean_pixel=mean(MassCenter_cell{1,i}{j}(1:pre_frames));
        plot(1:length(relevant_frames),MassCenter_cell{1,i}{j}-mean_pixel,'-')
        asparity_number=cat(1,asparity_number,[i,j]);
    end
end
title('\fontsize{16}Superposition of events - Mass center on X axis');
xlabel('Frames');
ylabel('normalized X_c [pixel]' );
ylim([-2.5 2.5])
legend(num2str(asparity_number))

%---plot2 - superposition of the movement in Y direction
figure
hold all
asparity_number=[];
for i=choose_videos
    for j=1:length(MassCenter_cell{2,i})
        mean_pixel=mean(MassCenter_cell{2,i}{j}(1:pre_frames));
        plot(1:length(relevant_frames),MassCenter_cell{2,i}{j}-mean_pixel,'-')
        asparity_number=cat(1,asparity_number,[i,j]);
    end
end
title('\fontsize{16}Superposition of events - Mass center on Y axis');
xlabel('Frames');
ylabel('normalized Y_c [pixel]' );
legend(num2str(asparity_number))

% 
% 
% %---plot3 - superposition of the speed in X direction
% figure
% hold all
% for i=1:3
% plot(1:length(relevant_frames)-1,MassCenter_cell{3,i},'-')
% end
% title('\fontsize{16}Superposition of events - speed of Mass center on X axis');
% xlabel('Frames');
% ylabel('speed [pixel/frames]' );
% 
% %---plot4 - superposition of the speed in Y direction
% figure
% hold all
% for i=1:3
% plot(1:length(relevant_frames)-1,MassCenter_cell{4,i},'-')
% end
% title('\fontsize{16}Superposition of events - speed of Mass center on Y axis');
% xlabel('Frames');
% ylabel('speed [pixel/frames]' );

% a=0;
% for i=1:number_of_videos
%     a=a+length(MassCenter_cell{3,i});
% end
% 
% A=zeros(1,a);
% for i=1:number_of_videos
%     for j=1:length(MassCenter_cell{3,i})
%         A(j+i-1)=max(MassCenter_cell{3,i}{j});
%     end
% end


% A=[];
% for i=1:number_of_videos
%     a=[];
%     for j=1:length(MassCenter_cell{3,i})
%         a(j)=max(MassCenter_cell{3,i}{j});
%     end
%     A=[A a];
% end