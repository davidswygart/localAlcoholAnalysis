function [spkCounts,  bpod] = binAroundTarget(c, target, binEdges, varargin)

spkCounts = nan(size(c,1), length(binEdges)-1);

for i=1:size(c,1)
    bpod = c.bpod(i);
    tVal = bpod.(target);

    %subtract off target value. There should be a better way of doing this.
    bpod.baselineStart = bpod.baselineStart - tVal;
    bpod.microInjectionStart = bpod.microInjectionStart - tVal;
    bpod.postInjectionStart = bpod.postInjectionStart - tVal;
    bpod.sipperStart = bpod.sipperStart - tVal;
    bpod.tailStart = bpod.tailStart - tVal;
    bpod.dropDeployed = bpod.dropDeployed - tVal;

    spk = c.spikeTimes{i};
    spk = spk-tVal;
    
    count = histcounts(spk,binEdges);

    if any(contains(varargin,'smooth'))
        meanISI = mean(diff(spk));
        binWidth = binEdges(2)-binEdges(1);
        count = imgaussfilt(count, meanISI/binWidth/4);
    end
    spkCounts(i,:) = count;
end
end
