

matFiles = string(ls('*.mat'));

%split files by experiment group
load('C:\Users\dis006\OneDrive - Indiana University\localAlcohol\group.mat')
controlFiles = matFiles(ismember(erase(matFiles, '.mat'), group.matName(contains(group.group,'control'))));
drinkFiles = matFiles(ismember(erase(matFiles, '.mat'), group.matName(contains(group.group,'drink'))));
injectedFiles = matFiles(ismember(erase(matFiles, '.mat'), group.matName(contains(group.group,'injected'))));

plotFiring(controlFiles)
plotFiring(drinkFiles)
plotFiring(injectedFiles)



function plotFiring(files)
binSize = 0.1; %bin width in seconds
maxTime = 47*60; %max time in seconds
xA = 0:binSize:maxTime;
nbins = length(xA)-1;
nFiles = size(files,1);

    for i=1:nFiles
        d = load(files(i,:));
    
        % choose clusters that pass quality metrics
        clusterInfo = d.clusterInfo;
        goodClusters = clusterInfo.fr>.1;
        goodClusters = goodClusters & ~contains(clusterInfo.phyLabel,'noise');
        goodClusters = goodClusters & clusterInfo.isiViolations<1;
        goodClusters = goodClusters & clusterInfo.presenceRatio>0.9;
        clusterInfo = clusterInfo(goodClusters,:);

        % bin spikes
    
    
        nclusters = size(clusterInfo,1);
        if nclusters<1
            break
        end


        bpodTimes = extractBpodTimes(d.events);
        spkCounts = nan(nclusters,nbins);
        plot(d.events.timestamp(d.events.line ==8), d.events.state(d.events.line ==8))
        for c=1:size(clusterInfo,1)
            spk = clusterInfo.spikeTimes{c};
            %spk = spk-bpodTimes.ba
            %spkCounts(c,:) = histcounts(spk, xA);
        end
    

    
    
    
    
    
    
    end

end
