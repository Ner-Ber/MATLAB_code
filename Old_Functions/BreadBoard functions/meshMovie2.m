function [Movie,Moviefigure] = meshMovie2(frames, Xmesh,Ymesh)
%frames should be a 3D matrix, 3rd dimension is time

Moviefigure=figure('Renderer','zbuffer');
mesh(Xmesh, Ymesh, frames(:,:,1));
axis tight
set(gca,'NextPlot','replaceChildren');
% Preallocate the struct array for the struct returned by getframe
Movie(size(frames,3)-0) = struct('cdata',[],'colormap',[]);
% Record the movie
for i = 1:size(frames,3)
    mesh(Xmesh, Ymesh, frames(:,:,i));
    Movie(i) = getframe;
end
end
% close