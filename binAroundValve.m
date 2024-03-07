function spkCounts = binAroundValve(c, binEdges,varargin)
ndrops = 60;
spkCounts = nan(size(c,1), length(binEdges)-1, ndrops);

for i=1:size(c,1)
    bpod = c.bpod(i);
    tVals = bpod.dropDeployed;

    for v=1:ndrops
        spk = c.spikeTimes{i} - tVals(v);
        if any(contains(varargin,'smooth'))
            binWidth = binEdges(2)-binEdges(1);
            npads = 5; %number of padding bins to put on each end, prevents end effects caused by smoothing
            padStart = binEdges(1)-binWidth*npads:binWidth:binEdges(1)-binWidth;
            padEnd =  binEdges(end)+binWidth:binWidth:binEdges(end)+binWidth*npads;
            newBins = [padStart, binEdges, padEnd];

            count = histcounts(spk,newBins);
            meanISI = mean(diff(spk));        
            count = imgaussfilt(count, meanISI/binWidth/4,'Padding','circular');
            count = count(npads+1:end-npads); % throw out padded bins
        else
            count = histcounts(spk,binEdges);
        end

        spkCounts(i,:,v) = count;
    end
end
end
