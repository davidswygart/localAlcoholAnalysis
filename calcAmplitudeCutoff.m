function clusters = calcAmplitudeCutoff(clusters)
    clusters.ampCutoff = cellfun(@cutoffFormula, clusters.allAmps);
end


function fraction_missing = cutoffFormula(amps)
    amps = abs(amps)*-1;
    nbins = 500;
    smoothing = 3;
    [h,b] = histcounts(amps,nbins);
    h = imgaussfilt(h,smoothing);
    pdf = h/ sum(h);
    [~, peakInd] = max(pdf);
    
    [~,minAfter] = min(abs(pdf(peakInd:end) - pdf(1)));
    G = minAfter+peakInd;
    binSize = mean(diff(b));
    
    fraction_missing = sum(pdf(G:end)) * binSize;
    fraction_missing = min([fraction_missing, 0.5]);
end