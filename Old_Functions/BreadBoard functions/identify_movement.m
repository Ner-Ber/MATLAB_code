%% parameters

%---generaye RGB...:
%# when no movie_name or movieType are inserted no movie will be generated, only frames.
movie_name = '';
movieType = '';
folderPath = 'C:\Users\owner\Documents\Shahar_Neri\friction_videos\front250mps_C001H001S0001';

%---crop frames
crop_movie = 1;
crop_region = 3:7;

%---crop relevant time
relevant_frames = 1:110;


%---measurment parameters
ImageData = get_cih_data(folderPath);
dt = 1/ImageData.RecordRate_fps;
pixel_resolution = 3e-6;        %pixel capture this much meters


%---analyze single row:
row_number = 2;

%---parabolas
plot_parabolas = 0;
plot_vertex = 1;

%% generate RGB intensity movie
[fixed_frames, fixed_frames_RGB, all_frames] =...
    generateFixedMovie(folderPath, movie_name, movieType);

%% crop frames of movies
if crop_movie
    all_frames = all_frames(crop_region,:,:);
    fixed_frames = fixed_frames(crop_region,:,:);
    fixed_frames_RGB = fixed_frames_RGB(crop_region,:,:);
end

%% analyze single row
single_row_over_time = fixed_frames(row_number,:,:);
single_row_over_time = permute(single_row_over_time, [3 2 1]);  %rotate into a 2D matrix, 1st dimension is time, 2nd is x axis

%---define time scale:
frameScale = 1:size(single_row_over_time,1);
timeScale = (frameScale-frameScale(1)+ImageData.StartFrame)./ImageData.RecordRate_fps;

frame_range = frameScale(relevant_frames);
time_range = timeScale(relevant_frames);

%---length scale
length_range = (1:size(single_row_over_time,2))*pixel_resolution;


num_of_frames = length(frame_range);
plot_colors = MyVaryColor(num_of_frames);

row_plot_in_time = figure;
x_axis = (1:size(single_row_over_time, 2))*pixel_resolution;
hold on;
for i = 1:length(frame_range)
    plot(x_axis,single_row_over_time(frame_range(i),:),'Color',plot_colors(i,:));
end
xlabel('locatiom [m]');
ylabel('intensity [0,1]');
caxis([time_range(1) time_range(end)]);
colorbar;
hold off;


%% follow parabola
%---insert coordinate ranges
Xranges = {...
    648:653;...
    280:286;...
    391:403;...
    615:623;...
    584:589;...
    537:545;...
    431:436;...
    };          %range of x coordinates. each cell indicates...
                %a different range which contains an asparity
[~, Xorder] = sort(cellfun(@(a) a(1), Xranges));
Xranges = Xranges(Xorder);


%---create cell to contain all parabola fits:
parabola_fits = cell(length(Xranges),length(frame_range));    %rows represent certain asperity, column is time step.
parb_vertex = zeros(size(parabola_fits));
figure(row_plot_in_time);
x_vertex = @(a,b) -b/(2*a);     % vertex x coor caculation
hold on;
for coor_set = 1:size(parabola_fits,1)
    for timeStep = 1:size(parabola_fits,2)
        %--- pick 3 highest points to fit the parabola
        [sorted, order] = sort(single_row_over_time(timeStep,Xranges{coor_set}));
        maximal_coordinate =Xranges{coor_set}(order(end));
        
        %---define number of dots to fit to only odd nums. even will vo
        %"rounded" up)
        number_of_fit_points = 3;
        x_for_fit = ((maximal_coordinate-floor(number_of_fit_points/2)):...
            (maximal_coordinate+floor(number_of_fit_points/2)));
        y_for_fit = single_row_over_time(timeStep,x_for_fit);
        
        %---fit selected points:
        p = polyfit(x_for_fit*pixel_resolution, y_for_fit, 2);   %fit a parabola
        
        %---plot parabola:  **THIS STEP IS EXTREMELY RESOURCE CONSUMING**
        if plot_parabolas
            x_parb = pixel_resolution*linspace(Xranges{coor_set}(1),Xranges{coor_set}(end),100);
            plot(x_parb,p(1).*x_parb.^2 ...
            +p(2).*x_parb ...
            +p(3),...
            '-.','Color', plot_colors(timeStep,:));
        end
        
        %---save the parabola:
        parabola_fits{coor_set,timeStep} = p;
        
        %---caclulate vertex 
        parb_vertex(coor_set,timeStep) = x_vertex(p(1),p(2));
        if plot_vertex
            line([parb_vertex(coor_set,timeStep), parb_vertex(coor_set,timeStep)], [0 0.6],...
                'Color',plot_colors(timeStep,:));
        end
        
    end
end


%---plot vertex speed
vertex_colors = MyVaryColor(length(Xranges));
for coor_set = 1:length(Xranges)                  %plot indicating dot on asparity plot
    plot(parb_vertex(coor_set, round(size(parb_vertex,2))), 0.6,...
        'o','MarkerFaceColor',vertex_colors(coor_set,:),'MarkerSize',10,'MarkerEdgeColor','none');
end

vertex_speed = diff(parb_vertex,1,2);   %in meters per frame
vertex_speed = vertex_speed*ImageData.RecordRate_fps;
vertexSpeedPlot = figure;
hold on;
for coor_set = 1:length(Xranges)
    plot(time_range(1:end-1), vertex_speed(coor_set,:),'Color',vertex_colors(coor_set,:));
end
hold off;

