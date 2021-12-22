function [h, hcb] = imagescWithNan(a,cm,nanclr,varargin)
% [h, hcb] = imagescWithNan(a,cm,nanclr,x,y)
%
% x,y are optional 1x2 vectors indicating the edges of the x axis and y
% axis. 
% IMAGESC with NaNs assigning a specific color to NaNs

[x,y,showNanInCbar] = setDefaults4function(varargin,[0,size(a,2)],[0,size(a,1)],0);

%# find minimum and maximum
amin=min(a(:));
amax=max(a(:));
%# size of colormap
n = size(cm,1);
%# color step
dmap=(amax-amin)/n;

%# standard imagesc
him = imagesc(x,y,a);
%# add nan color to colormap
colormap([nanclr; cm]);
%# changing color limits
if ~showNanInCbar
    caxis([amin-dmap amax]);
    %# place a colorbar
    hcb = colorbar;
else
    hcb = colorbar;
end
%# change Y limit for colorbar to avoid showing NaN color
ylim(hcb,[amin amax])

if nargout > 0
    h = him;
end

end