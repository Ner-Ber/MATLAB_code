function CSV_fig = plot_from_CSV(CSV_path)


[num,txt,raw] = xlsread(CSV_path);
data_starting_row = 18;

measure_data = num(data_starting_row:end,:);

CSV_fig = figure;
hold all;
plot(measure_data(:,1),measure_data(:,2),'.');
plot(measure_data(:,1),measure_data(:,3),'.');
plot(measure_data(:,1),measure_data(:,4),'.');



end