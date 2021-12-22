start=35;  %%% the exact location of the event   
stop=50;
line_number=1;  %%%choose line of pixels (out of 8) to plot
frames=fixed_frames; %%%choose fixed_frames / all_frames (with blinking)
figure
hold all
for i=start:1:stop   %% choose jumps of 2 to get rid of the blinking
    x=1:size(all_frames,2);
     plot(x,frames(line_number,:,i))
end
hold off
axis([340 390 0 0.6])