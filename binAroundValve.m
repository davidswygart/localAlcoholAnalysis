function spkCounts = binAroundValve(c, binEdges,varargin)
ndrops = 60;
spkCounts = nan(size(c,1), length(binEdges)-1, ndrops);

for i=1:size(c,1)
    bpod = c.bpod(i);
    tVals = bpod.dropDeployed;

    for v=1:ndrops
        spk = c.spikeTimes{i} - tVals(v);
        count = histcounts(spk,binEdges);
        if any(contains(varargin,'smooth'))
            meanISI = mean(diff(spk));
            binWidth = binEdges(2)-binEdges(1);
            count = imgaussfilt(count, meanISI/binWidth/4);
        end
        spkCounts(i,:,v) = count;
    end
end
end
