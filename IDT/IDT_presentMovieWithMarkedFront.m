function IDT_presentMovieWithMarkedFront(exp_dir,eventNum)

[Mov4D_colorMat] = IDT_createNormalizedMovie(exp_dir,eventNum);
[uniqTimes, uniqLocs] = IDT_get2D_front(exp_dir,eventNum);


end