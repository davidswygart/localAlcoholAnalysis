%% ISI violations
figure(10);clf
binEdges = logspace(log10(1), log10(max(allClusters.isiViolations)+1),40);
histogram(allClusters.isiViolations+1, binEdges)
set(gca,'xscale','log')
hold on
plot([2,2], ylim, 'r--')
hold off
xlabel('ISI violations (Hill equation) + 1')
ylabel('Number of clusters')
oldX = xlim;
xlim([oldX(1), 2e3])
saveas(gcf,'isiViolations.svg')

%% Presence ratio
figure(11);clf
histogram(allClusters.presenceRatio,40)
hold on
plot([.9,.9], ylim, 'r--')
hold off
xlabel('Presence ratio')
ylabel('Number of clusters')
saveas(gcf, 'presenceRatio.svg')




%% survival N
figure(1)
plotSurvival(allClusters(contains(allClusters.group, 'control'),:))
ylim([0,1122])

saveas(gcf, 'clusterN_control.svg')

figure(2)
plotSurvival(allClusters(contains(allClusters.group, 'drink'),:))
ylim([0,1007])
saveas(gcf, 'clusterN_drink.svg')

figure(3)
plotSurvival(allClusters(contains(allClusters.group, 'inject'),:))
ylim([0,681])
saveas(gcf, 'clusterN_inject.svg')


function plotSurvival(clusters)
[uNames,~,ui] = unique(clusters.matName);
nclusters = nan(4,length(uNames));

for i=1:length(uNames)
    nclusters(1,i) = sum(ui==i);
    nclusters(2,i) = sum(ui==i & ~contains(clusters.phyLabel,'noise'));
    nclusters(3,i) = sum(ui==i & ~contains(clusters.phyLabel,'noise') & clusters.isiViolations<1);
    nclusters(4,i) = sum(ui==i & ~contains(clusters.phyLabel,'noise') & clusters.isiViolations<1 & clusters.presenceRatio>.9);
end


bar(nclusters, 'stacked')
ylabel('number of clusters')
set(gca,'xtick',[])
xlim([0.5,4.5])
end