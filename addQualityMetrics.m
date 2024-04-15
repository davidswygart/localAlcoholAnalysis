function clusters = addQualityMetrics(clusters)
    clusters.spikeTimes = cellfun(@sort, clusters.spikeTimes, 'UniformOutput', false);% sort spike times for each cluster
    clusters = calcIsiViolations(clusters);
    clusters = calcPresenceRatio(clusters);
    clusters = calcAmplitudeCutoff(clusters);
end