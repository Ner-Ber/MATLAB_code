start=1;  %%% the exact location of the event   
stop=30;
figure
hold all
for i=start:1:stop   %% jumps of 2 to get rid of the blinking
    x=1:size(RowOverTime,2);
    plot(x,RowOverTime(i,:))
end
hold off
