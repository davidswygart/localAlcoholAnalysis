
function loc3d = calc3Dposition(loc4d, angle)
loc3d = nan(1,3);

loc3d(1) = loc4d(1) + cosd(angle)*loc4d(4);
loc3d(2) = loc4d(2); % Y is the same in 3D as 4D
loc3d(3) = loc4d(3) - sind(angle)*loc4d(4);
end