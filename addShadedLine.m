function h = addShadedLine(x, ymat, color, name)
    ymat = squeeze(ymat);
    avg = mean(ymat,1);
    stdev = std(ymat,0, 1);
    n = size(ymat,1);
    sem = stdev / sqrt(n);
    if isempty(x)
        x = 1:length(avg);
    end
    h = shadedErrorBar(x, avg,sem, 'lineProps',color);
    h.mainLine.DisplayName = name;
end