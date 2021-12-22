function [Place,Angle,Length] = Movie_calcPlaceAngleLength()
% when marking pairs of data cursors multiple slopes (just starting and
% ending) this function will retrive the Place,Angle and Length of the
% slopes. 


d = data_tip_get_from_mesh;
xx = d.x;
yy = d.y;
[XX, order] = sort(xx);
YY = yy(order);

[~, I] = max([YY(1),YY(2)]);

if I==1
    XX_uppers = XX(1:2:end);
    XX_lowers = XX(2:2:end);
    YY_uppers = YY(1:2:end);
    YY_lowers = YY(2:2:end);
elseif I==2
    XX_uppers = XX(2:2:end);
    XX_lowers = XX(1:2:end);
    YY_uppers = YY(2:2:end);
    YY_lowers = YY(1:2:end);
end
Place = XX_lowers;
Angle = atand((XX_uppers-XX_lowers)./(YY_uppers-YY_lowers));
Length = sqrt((XX_uppers-XX_lowers).^2+(YY_uppers-YY_lowers).^2);




end