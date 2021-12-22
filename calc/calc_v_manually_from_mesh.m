function [v]=calc_v_manually_from_mesh
%The function gets the data cursors from the mesh and calculates the average velocity betwean the points. 

dp=data_tip_get_from_mesh;

[x,s]=sort(dp.x);
y=dp.y(s);
v=diff(x)./diff(y);
