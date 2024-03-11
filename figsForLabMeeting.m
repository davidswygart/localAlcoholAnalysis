%% filter out bad clusters
goodClusters = allClusters;
goodClusters = goodClusters(goodClusters.fr>0.1, :);
goodClusters = goodClusters(~contains(goodClusters.phyLabel,'noise'), :);
goodClusters = goodClusters(goodClusters.presenceRatio>0.9, :);
goodClusters = goodClusters(goodClusters.isiViolations<1, :);

%% Give example heatmaps for specific mice (spike count)
binWidth=0.1;
binEdges = -10*60:binWidth:60*37;
target = 'microInjectionStart';
cLabel = 'spikes';
cRange = [0,8];

figure(1); clf
t = tiledlayout(3,1,'TileSpacing','Compact','Padding','Compact');

nexttile
exampleControl = goodClusters(contains(goodClusters.matName, '2023-11-16_10-56-52'), :);
[spkCounts,  bpod]  = binAroundTarget(exampleControl, target, binEdges);
plotSpikeHeatmapWithBpod(spkCounts, binEdges, bpod, cLabel,cRange)
set(gca,'xticklabel',{[]})
title('Control')

nexttile
exampleDrink = goodClusters(contains(goodClusters.matName, '2023-11-15_12-27-09'), :);
[spkCounts,  bpod]  = binAroundTarget(exampleDrink, target, binEdges);
plotSpikeHeatmapWithBpod(spkCounts, binEdges, bpod, cLabel,cRange)
set(gca,'xticklabel',{[]})
title('Drink')

nexttile
exampleInj = goodClusters(contains(goodClusters.matName, '2023-11-16_13-16-55'), :);
[spkCounts,  bpod]  = binAroundTarget(exampleInj, target, binEdges);
plotSpikeHeatmapWithBpod(spkCounts, binEdges, bpod, cLabel,cRange)
xlabel('Time (0.1s bins)')
title('Inject')

saveas(gcf, 'example_spikes.svg')

%% Give example heatmaps for specific mice (spike count - smoothed)
binWidth=0.1;
binEdges = -10*60:binWidth:60*37;
target = 'microInjectionStart';
cLabel = 'smoothed spikes';
cRange = [0,8];

figure(2); clf
t = tiledlayout(3,1,'TileSpacing','Compact','Padding','Compact');

nexttile
exampleControl = goodClusters(contains(goodClusters.matName, '2023-11-16_10-56-52'), :);
[spkCounts,  bpod]  = binAroundTarget(exampleControl, target, binEdges, 'smooth');
plotSpikeHeatmapWithBpod(spkCounts, binEdges, bpod, cLabel,cRange)
set(gca,'xticklabel',{[]})
title('Control')

nexttile
exampleDrink = goodClusters(contains(goodClusters.matName, '2023-11-15_12-27-09'), :);
[spkCounts,  bpod]  = binAroundTarget(exampleDrink, target, binEdges, 'smooth');
plotSpikeHeatmapWithBpod(spkCounts, binEdges, bpod, cLabel,cRange)
set(gca,'xticklabel',{[]})
title('Drink')

nexttile
exampleInj = goodClusters(contains(goodClusters.matName, '2023-11-16_13-16-55'), :);
[spkCounts,  bpod]  = binAroundTarget(exampleInj, target, binEdges, 'smooth');
plotSpikeHeatmapWithBpod(spkCounts, binEdges, bpod, cLabel,cRange)
xlabel('Time (0.1s bins)')
title('Inject')

saveas(gcf, 'example_spikes_smoothed.svg')

%% Give example heatmaps for specific mice (zscore)
binWidth=0.1;
binEdges = -10*60:binWidth:60*37;
target = 'microInjectionStart';
cLabel = 'zscore';
cRange = [-.5,8];

figure(3); clf
t = tiledlayout(3,1,'TileSpacing','Compact','Padding','Compact');

nexttile
exampleControl = goodClusters(contains(goodClusters.matName, '2023-11-16_10-56-52'), :);
[spkCounts,  bpod]  = binAroundTarget(exampleControl, target, binEdges, 'smooth');
spkCounts = zscore(spkCounts,0,2);
plotSpikeHeatmapWithBpod(spkCounts, binEdges, bpod, cLabel,cRange)
set(gca,'xticklabel',{[]})
title('Control')

nexttile
exampleDrink = goodClusters(contains(goodClusters.matName, '2023-11-15_12-27-09'), :);
[spkCounts,  bpod]  = binAroundTarget(exampleDrink, target, binEdges, 'smooth');
spkCounts = zscore(spkCounts,0,2);
plotSpikeHeatmapWithBpod(spkCounts, binEdges, bpod, cLabel,cRange)
set(gca,'xticklabel',{[]})
title('Drink')

nexttile
exampleInj = goodClusters(contains(goodClusters.matName, '2023-11-16_13-16-55'), :);
[spkCounts,  bpod]  = binAroundTarget(exampleInj, target, binEdges, 'smooth');
spkCounts = zscore(spkCounts,0,2);
plotSpikeHeatmapWithBpod(spkCounts, binEdges, bpod, cLabel,cRange)
xlabel('Time (0.1s bins)')
title('Inject')

saveas(gcf, 'example_zscore.svg')

%% spikes around injection 
binWidth=10;
binEdges = -200:binWidth:60*12;
target = 'microInjectionStart';
cLabel = 'zscore';
cRange = [-.5,8];

[spkCounts,  bpod]  = binAroundTarget(goodClusters, target, binEdges, 'smooth');
spkZ = zscore(spkCounts,0, 2);

isControl = contains(goodClusters.group,'control') | contains(goodClusters.group,'drink');
control = spkZ(isControl,:);
isInject = contains(goodClusters.group,'inject');
inject = spkZ(isInject,:);

%% spikes around injection (heatmap)
figure(4); clf
t = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');

nexttile
plotSpikeHeatmapWithBpod(control, binEdges, bpod, cLabel,cRange)
set(gca,'xticklabel',{[]})
title('Control')

nexttile
plotSpikeHeatmapWithBpod(inject, binEdges, bpod, cLabel,cRange)
title('Inject')
xlabel('Time (10s bins)')


%% spikes around injection (average)
figure(5); clf
x = binEdges(1:end-1);
addShadedLine(x,control,'b','Control')
hold on
addShadedLine(x,inject,'r','Inject')

plotBpod(bpod)
xlim([binEdges(1), binEdges(end)])

h = multipleTtest({control,inject});

y = ylim();
y = ones(length(x),1) * y(1);
hold on
scatter(x(h),y(h),'*k')
legend('Control','Inject')
xlabel('time (s)')
ylabel('zscore')

%% spikes around injection (histogram of specific points)
figure(6); clf
t = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
nbins = 50;

nexttile
t = 40;
[~,tInd] = min(abs(binEdges - t));
histogram(control(:,tInd),nbins)
hold on
histogram(inject(:,tInd),nbins)
legend('control','inject')
ylabel('number of clusters')
title(['Timepoint: ', num2str(t), 's'])

nexttile
t = 710;
[~,tInd] = min(abs(binEdges - t));
histogram(control(:,tInd),nbins)
hold on
histogram(inject(:,tInd),nbins)
legend('control','inject')
ylabel('number of clusters')
title(['Timepoint: ', num2str(t), 's'])
xlabel('zscore')


