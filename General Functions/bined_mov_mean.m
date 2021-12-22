function [x_out, mean_out] = bined_mov_mean(x,v,bin_w)
    
    bin_edges = (min(x)-bin_w/2):bin_w:(max(x)+bin_w/2);
    Y = discretize(x,bin_edges);
    
end