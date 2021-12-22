function [dp]=data_tip_get_from_mesh_all_fields
%The function gets the data cursors from the mesh and calculates the average velocity betwean the points. 

dcm_obj = datacursormode(gcf);
info_struct = getCursorInfo(dcm_obj);

number_of_cursors=length(info_struct);
for j=1:number_of_cursors
dp.x(j)=info_struct(j).Position(1);
dp.y(j)=info_struct(j).Position(2);
% dp.z(j)=info_struct(j).Position(3);
end
