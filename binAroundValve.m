function spkCounts = binAroundValve(c, binEdges)
ndrops = 60;
spkCounts = nan(size(c,1), length(binEdges)-1, ndrops);

for i=1:size(c,1)
    bpod = c.bpod(i);
    tVals = bpod.dropDeployed;

    for v=1:ndrops
        spk = c.spikeTimes{i} - tVals(v);
        spkCounts(i,:,v) = histcounts(spk,binEdges);
    end
end
end
