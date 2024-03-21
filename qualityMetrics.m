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

%% Depth cuttoff
c = allClusters(~contains(allClusters.phyLabel,'noise') & allClusters.isiViolations<1 & allClusters.presenceRatio>.9, :);
figure(1);clf
scatter(c.depth,c.distFromInj, 'b.')
ylim([0,3000])
xlabel('Distance from NP shank tip (um)')
ylabel('Distance from injection (um)')
hold on
plot([1100,1100], ylim, 'k--')

%% survival N
figure(1); clf
plotSurvival(allClusters(contains(allClusters.group, 'control'),:))
saveas(gcf, 'clusterN_control.svg')

figure(2); clf
plotSurvival(allClusters(contains(allClusters.group, 'drink'),:))
saveas(gcf, 'clusterN_drink.svg')

figure(3); clf
plotSurvival(allClusters(contains(allClusters.group, 'inject'),:))
saveas(gcf, 'clusterN_inject.svg')

%% Equal depth?
injDist = goodClusters.depth(contains(goodClusters.group,'inject'));
controlDist = goodClusters.depth(contains(goodClusters.group,'control') | contains(goodClusters.group,'drink'));


figure(4); clf
[f,x] = ecdf(controlDist);
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',2,'Color','b')
hold on
[f,x] = ecdf(injDist);
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',2,'Color','r')

plot([0,0],[0,0.5], '--k')
xlabel('Distance from tip of probe')
ylabel('Folded probability')
legend("Control", "Inject","Location","northeast")

%xlim([600,1500])

[h,p] = kstest2(controlDist,injDist);
text(820,.25,['p=',num2str(p,2)])

%% Equal distance to injection site?
injDist = goodClusters.distFromInj(contains(goodClusters.group,'inject'));
controlDist = goodClusters.distFromInj(contains(goodClusters.group,'control') | contains(goodClusters.group,'drink'));


figure(4); clf
[f,x] = ecdf(controlDist);
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',2,'Color','b')
hold on
[f,x] = ecdf(injDist);
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',2,'Color','r')

plot([0,0],[0,0.5], '--k')
xlabel('Distance from injection')
ylabel('Folded probability')
legend("Control", "Inject","Location","northwest")

xlim([600,1500])

[h,p] = kstest2(controlDist,injDist);
text(620,.25,['p=',num2str(p,2)])

%%
function plotSurvival(clusters)
[uNames,~,ui] = unique(clusters.matName);
nclusters = nan(5,length(uNames));

for i=1:length(uNames)
    nclusters(1,i) = sum(ui==i);
    nclusters(2,i) = sum(ui==i & ~contains(clusters.phyLabel,'noise'));
    nclusters(3,i) = sum(ui==i & ~contains(clusters.phyLabel,'noise') & clusters.isiViolations<1);
    nclusters(4,i) = sum(ui==i & ~contains(clusters.phyLabel,'noise') & clusters.isiViolations<1 & clusters.presenceRatio>.9);
    nclusters(5,i) = sum(ui==i & ~contains(clusters.phyLabel,'noise') & clusters.isiViolations<1 & clusters.presenceRatio>.9 & clusters.depth < 1100);
end


bar(["FR>0.1Hz","Curation","ISI violation","Presence ratio","Depth<1100"],nclusters, 'stacked')
ylabel('number of clusters')
%set(gca,'xtick',[])
%xlim([0.5,4.5])
ylim([0,sum(nclusters(1,:))])

c = clusters(~contains(clusters.phyLabel,'noise') & clusters.isiViolations<1 & clusters.presenceRatio>.9 & clusters.depth < 1100, :);
nMouse = length(unique(c.matName));
nCluster_total = size(c,1);

msg = sprintf("Clusters=%i",nCluster_total);
msg = [msg; sprintf("Mice=%i",nMouse)];
text(5,nCluster_total*4,msg,"HorizontalAlignment","center")

end


