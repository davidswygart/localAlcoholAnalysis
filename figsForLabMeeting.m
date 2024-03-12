%% filter out bad clusters
% goodClusters = allClusters;
% goodClusters = goodClusters(goodClusters.fr>0.1, :);
% goodClusters = goodClusters(~contains(goodClusters.phyLabel,'noise'), :);
% goodClusters = goodClusters(goodClusters.presenceRatio>0.9, :);
% goodClusters = goodClusters(goodClusters.isiViolations<1, :);
% goodClusters = sortrows(goodClusters, "fr", 'descend');

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

saveas(gcf, 'injection_heatmap.svg')
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

xlabel('time (s)')
ylabel('zscore')
plot(xlim(), [0,0], '--k')

legend('Control','Inject')

saveas(gcf, 'injection_mean.svg')
%% spikes around injection (histogram of specific points)
figure(6); clf
t = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
binEdges_pntCompare = -2.4 :0.3: 5.4;

nexttile
t = 40;
[~,tInd] = min(abs(binEdges - t));
histogram(control(:,tInd),binEdges_pntCompare)
hold on
histogram(inject(:,tInd),binEdges_pntCompare)
legend('control','inject')
ylabel('number of clusters')
title(['Timepoint: ', num2str(t), 's'])

nexttile
t = 710;
[~,tInd] = min(abs(binEdges - t));
histogram(control(:,tInd),binEdges_pntCompare)
hold on
histogram(inject(:,tInd),binEdges_pntCompare)
legend('control','inject')
ylabel('number of clusters')
title(['Timepoint: ', num2str(t), 's'])
xlabel('zscore')

saveas(gcf, 'injection_timepointHistogram.svg')

%% Spikes Around Drinking - heatmap
binWidth=0.1;
binEdges = -200:binWidth:60*22;
target = 'sipperStart';
cLabel = 'zscore';
cRange = [-.5,8];

[spkCounts,  bpod]  = binAroundTarget(goodClusters, target, binEdges, 'smooth');
spkZ = zscore(spkCounts,0, 2);

isControl = contains(goodClusters.group,'control');
control = spkZ(isControl,:);
isDrink = contains(goodClusters.group,'drink');
drink = spkZ(isDrink,:);
isInject = contains(goodClusters.group,'inject');
inject = spkZ(isInject,:);

figure(8); clf
t = tiledlayout(3,1,'TileSpacing','Compact','Padding','Compact');

nexttile
plotSpikeHeatmapWithBpod(control, binEdges, bpod, cLabel,cRange)
set(gca,'xticklabel',{[]})
title('Control')

nexttile
plotSpikeHeatmapWithBpod(drink, binEdges, bpod, cLabel,cRange)
title('Drink')

nexttile
plotSpikeHeatmapWithBpod(inject, binEdges, bpod, cLabel,cRange)
title('Inject')
xlabel('Time (0.1s bins)')

saveas(gcf, 'Sipping_heatmap.svg')

%% Spikes around Drinking - mean
figure(9); clf
x = binEdges(1:end-1);
addShadedLine(x,control,'b','Control')
hold on
addShadedLine(x,drink,'r','Drink')
addShadedLine(x,inject,'g','Inject')

plotBpod(bpod)
xlim([binEdges(1), binEdges(end)])

legend('Control','Drink','Inject')
xlabel('time (s)')
ylabel('zscore')

saveas(gcf, 'Sipper_mean.svg')

%% Spikes around sipper valve
binWidth=0.1;
binEdges = -1:binWidth:14;
spkCounts  = binAroundValve(goodClusters, binEdges,'smooth');
spkAvg = mean(spkCounts ,3);
spkZ = zscore(spkAvg,0, 2);

isControl = contains(goodClusters.group,'control');
control = spkZ(isControl,:);
isDrink = contains(goodClusters.group,'drink');
drink = spkZ(isDrink,:);
isInject = contains(goodClusters.group,'inject');
inject = spkZ(isInject,:);

%% Spikes around sipper valve - heatmap
figure(10); clf
t = tiledlayout(3,1,'TileSpacing','Compact','Padding','Compact');

nexttile
plotSpikeHeatmapWithBpod(control, binEdges, [], cLabel,cRange)
x = [(0-binEdges(1))/binWidth, (10.1-binEdges(1))/binWidth];
scatter(x, [0,0], 'r*')
yL = ylim();
ylim([0,yL(2)])
set(gca,'xticklabel',{[]})
title('Control')

nexttile
plotSpikeHeatmapWithBpod(drink, binEdges, [], cLabel,cRange)
scatter(x, [0,0], 'r*')
yL = ylim();
ylim([0,yL(2)])
set(gca,'xticklabel',{[]})
title('Drink')

nexttile
plotSpikeHeatmapWithBpod(inject, binEdges, [], cLabel,cRange)
scatter(x, [0,0], 'r*')
yL = ylim();
ylim([0,yL(2)])
title('Inject')
xlabel('Time (0.1s bins)')

saveas(gcf, 'valvesHeatmap.svg')

%% Spikes around sipper valve - mean
figure(11); clf
title('Average spiking after valve open')
hold on
x = binEdges(1:end-1);
addShadedLine(x, control, 'b','Control')
addShadedLine(x, drink, 'r', 'Drink')
addShadedLine(x, inject, 'g', 'Inject')
plot(xlim(), [0,0], '--k')

legend('Control', 'Drink', 'Inject')
ylabel('zscore')
xlabel('time (s)')
saveas(gcf, 'valvesMean.svg')
