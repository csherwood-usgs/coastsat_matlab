function [x, y] = transect2utm(x_on,y_on,Y,angle)
% [x, y] = transect2utm(x_on,y_on,Y,angle)
% Convert an array of transect Y values to UTM
x=x_on+Y.*cosd(angle);
y=y_on+Y.*sind(angle);
return