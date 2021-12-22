function [] = Movie_stack_fit_parabolas(dataMatrix)


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