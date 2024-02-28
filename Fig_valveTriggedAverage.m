%load("allClusters.mat")

allClusters = allClusters(allClusters.fr>0.1, :);
allClusters = allClusters(~contains(allClusters.phyLabel,'noise'), :);
allClusters = allClusters(allClusters.presenceRatio>0.9, :);
allClusters = allClusters(allClusters.isiViolations<1, :);

%% choose bin width
binWidth=0.1;
binEdges = -2:binWidth:15;

%%
includedValves = 20:40;
spkCounts  = binAroundValve(allClusters, binEdges);
spkCounts = mean(spkCounts(:,:,includedValves) ,3);

spkRate = spkCounts / binWidth;
spkZ = zscore(spkCounts,0, 2);

yname = 'zscore';
control = spkZ(contains(allClusters.group,'control'), :);
drink  = spkZ(contains(allClusters.group,'drink'), :);
inject = spkZ(contains(allClusters.group,'injected'), :);


%%
avgControl = mean(control,1);
stdControl = std(control,0, 1);
semControl = stdControl / sqrt(size(control,1));

avgDrink = mean(drink,1);
stdDrink = std(drink,0, 1);
semDrink = stdDrink / sqrt(size(drink,1));

avgInject = mean(inject,1);
stdInject = std(inject,0, 1);
semInject = stdInject / sqrt(size(inject,1));

figure(4)
clf
shadedErrorBar(binEdges(1:end-1), avgControl,semControl, 'lineProps', 'b')
hold on
shadedErrorBar(binEdges(1:end-1), avgDrink,semDrink, 'lineProps', 'r')
shadedErrorBar(binEdges(1:end-1), avgInject,semInject, 'lineProps', 'g')

legend('control','drink','inject')
ylabel(yname)
xlabel('time (s)')
xlim([binEdges(1), binEdges(end)])
hold off