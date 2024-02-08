function clusters = calcPresenceRatio(clusters)

%% bin spikes
binSize = 60; %bin width in seconds
maxTime = max(cellfun(@max, clusters.spikeTimes)); %max time in seconds
binEdges = 0:binSize:maxTime;
nbins = length(binEdges)-1;
nclusters = size(clusters,1);
spkCounts = nan(nclusters,nbins);
for i =1:nclusters
    spkCounts(i,:) = histcounts(clusters.spikeTimes{i}, binEdges);
end

%% find proportio of binds with at least 1 spike
clusters.presenceRatio = mean(spkCounts>0,2);


end

