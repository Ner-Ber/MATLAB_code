function Mov = meshMovie(frames, Xmesh,Ymesh)
%frames should be a 3D matrix, 3rd dimension is time

figure;
Mov = [];
for i = 1:size(frames,3)
    mesh(Xmesh, Ymesh, frames(:,:,i));
    F = getframe;
    try
        Mov = cat(4, Mov, F.cdata);
    catch
        return;
    end
end