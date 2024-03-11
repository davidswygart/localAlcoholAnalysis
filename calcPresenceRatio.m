function clusters = calcPresenceRatio(clusters)
 [~, ~, datasetInds]= unique(clusters.matName);
 clusters.presenceRatio(:) = nan;

for i = 1:max(datasetInds)
    currentDataset = datasetInds == i;
    binSize = 60; %bin width in seconds
    spkTimes = clusters.spikeTimes(currentDataset);
    minTime = min(cellfun(@min, spkTimes)); %min time for this dataset
    maxTime = max(cellfun(@max, spkTimes)); %max time for this dataset
    binEdges = minTime:binSize:maxTime;

    p = cellfun(@(x) presence(x, binEdges), spkTimes);
    clusters.presenceRatio(currentDataset) = p;
end
end

function p = presence(spkTimes,binEdges)
    counts = histcounts(spkTimes, binEdges);
    p = mean(counts>0);
end

