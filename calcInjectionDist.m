function clusters = calcInjectionDist(clusters)
angle = 75;

needleOD = 210; % 33G OD in microns
needleBevel = 12; % standard bevel angle of Hamilton Point style 4
needle_yOffset = needleOD/2;
needle_dOffset = needleOD / tand(needleBevel) / 2;

for i = 1:size(clusters,1)
    needle4D = clusters.needleLoc4D{i};
    needle4D(2) = needle4D(2) + needle_yOffset;
    needle4D(4) = needle4D(4) - needle_dOffset;
    needle3D = calc3Dposition(needle4D, angle);

    np4D = clusters.npLoc4D{i};
    np4D(4) = np4D(4) - clusters.depth(i);
    np3D = calc3Dposition(np4D, angle);
    np3D(1) = np3D(1)*-1;
    
    clusters.distFromInj(i) = sqrt(sum((needle3D - np3D).^2));
end
end