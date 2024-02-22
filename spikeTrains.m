matFiles = string(ls('*.mat'));

%split files by experiment group
load('C:\Users\david\OneDrive - Indiana University\localAlcohol\ephysData\group.mat')
controlFiles = matFiles(ismember(erase(matFiles, '.mat'), group.matName(contains(group.group,'control'))));
drinkFiles = matFiles(ismember(erase(matFiles, '.mat'), group.matName(contains(group.group,'drink'))));
injectedFiles = matFiles(ismember(erase(matFiles, '.mat'), group.matName(contains(group.group,'injected'))));

figure(1)
[spk, bpod] = getbins(controlFiles);
plotTrainsOneHeatmap(spk,bpod)
title('Control')

figure(2)
[spk, bpod] = getbins(drinkFiles);
plotTrainsOneHeatmap(spk,bpod)
title('Drink alcohol')

figure(3)
[spk, bpod] = getbins(injectedFiles);
plotTrainsOneHeatmap(spk,bpod)
title('inject alcohol')

function [allSpkCounts, allbpodTimes] = getbins(files)

nFiles = size(files,1);

allSpkCounts = {};
allbpodTimes = {};
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


        binSize = 0.1; %bin width in seconds
        startTime = 0;%-5*60; 
        endTime = 47*60; %max time in seconds
        binEdges = startTime:binSize:endTime;
        nbins = length(binEdges)-1;
        spkCounts = nan(nclusters,nbins);
        

        bpodTimes = extractBpodTimes(d.events);
        for c=1:size(clusterInfo,1)
            spk = clusterInfo.spikeTimes{c};
            spk = spk - bpodTimes.baselineStart;
            spkCounts(c,:) = histcounts(spk,binEdges);
        end
        imagesc(spkCounts)
        allSpkCounts = [allSpkCounts; spkCounts];

        
        bpodTimes.dropDeployed = bpodTimes.dropDeployed - bpodTimes.baselineStart;
        bpodTimes.microInjectionStart = bpodTimes.microInjectionStart - bpodTimes.baselineStart;
        bpodTimes.postInjectionStart = bpodTimes.postInjectionStart - bpodTimes.baselineStart;
        bpodTimes.sipperStart = bpodTimes.sipperStart - bpodTimes.baselineStart;
        bpodTimes.tailStart = bpodTimes.tailStart - bpodTimes.baselineStart;
        bpodTimes.baselineStart = bpodTimes.baselineStart - bpodTimes.baselineStart;
        allbpodTimes = [allbpodTimes; bpodTimes];
    end
    %imagesc(allSpkCounts)
end

function plotTrain(spk, bpod)
clf
nDatasets = size(spk,1);
tiledlayout(nDatasets,1,"TileSpacing","compact")

for i=1:nDatasets
    %subplot(nDatasets,1,i)
    nexttile
    imagesc(spk{i})
    set(gca,'xticklabel',{[]})
end
end


function plotTrainsOneHeatmap(spk, bpod)

nbins = size(spk{1},2);
nNanRows = 2;
spkMat = cellfun(@(x) [nan(nNanRows,nbins);x],spk,'UniformOutput',false);
spkMat = cell2mat(spkMat);
spk

clf
h = imagesc(spkMat);
set(h, 'AlphaData', ~isnan(spkMat))

head = find(isnan(spkMat(:,1)));
head = head(1:nNanRows:end);

hold on
for i=1:length(head)
    y = [head(i),head(i)+1];
    b = bpod{i};

    x = b.microInjectionStart*10;
    plot([x,x], y,'r')
    text(x,y(1), 'inj')
    x = b.sipperStart*10;
    plot([x,x], y, 'r')
    text(x,y(1), 'sip')
    x = b.tailStart*10;
    plot([x,x], y, 'r')
    text(x,y(1), 'tail')
    x = b.postInjectionStart*10;
    plot([x,x], y, 'r')
    text(x,y(1), 'post-inj')

    x = b.dropDeployed*10;
    y = ones(length(x),1) * head(i)+1; 
    scatter(x,y, 'r*')
    
end
hold off
xlabel('time (s*10)')
ylabel('cluster')
c = colorbar;
ylabel(c,'spike counts')
end

