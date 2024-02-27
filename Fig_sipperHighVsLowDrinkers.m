%load("allClusters.mat")

allClusters = allClusters(allClusters.fr>0.1, :);
allClusters = allClusters(~contains(allClusters.phyLabel,'noise'), :);
allClusters = allClusters(allClusters.presenceRatio>0.9, :);
allClusters = allClusters(allClusters.isiViolations<1, :);

%% choose anchor point and bin width
binWidth=10;

binEdges = -200:binWidth:60*25;
target = 'sipperStart';

% binEdges = -600:binWidth:60*12;
% target = 'microInjectionStart';
%%

[spkCounts,  bpod]  = binAroundTarget(allClusters, target, binEdges);
spkRate = spkCounts / binWidth;
spkZ = zscore(spkCounts,0, 2);

 %% Choose spike rate or z score
yname = 'zscore';
control = spkZ(contains(allClusters.group,'control'), :);

density = 1.006;
concentrationAlcohol = 0.2;
threshold = 2  /concentrationAlcohol/density % ml/kg threshold
drinkHigh = spkZ(contains(allClusters.group, 'drink') & allClusters.mlPerkg > threshold,:);
drinkLow = spkZ(contains(allClusters.group, 'drink') & allClusters.mlPerkg < threshold,:);

% threshold = 20 % mg% threshold
% drinkHigh = spkZ(contains(allClusters.group, 'drink') & allClusters.brainAlcohol > threshold,:);
% drinkLow = spkZ(contains(allClusters.group, 'drink') & allClusters.brainAlcohol < threshold,:);
%%
avgControl = mean(control,1);
stdControl = std(control,0, 1);
semControl = stdControl / sqrt(size(control,1));

avgHigh = mean(drinkHigh,1);
stdHigh = std(drinkHigh,0, 1);
semHigh = stdHigh / sqrt(size(drinkHigh,1));

avgLow = mean(drinkLow,1);
stdLow = std(drinkLow,0, 1);
semLow = stdLow / sqrt(size(drinkLow,1));

figure(3)
clf
shadedErrorBar(binEdges(1:end-1), avgControl,semControl, 'lineProps', 'b')
hold on
shadedErrorBar(binEdges(1:end-1), avgHigh,semHigh, 'lineProps', 'r')
shadedErrorBar(binEdges(1:end-1), avgLow,semLow, 'lineProps', 'g')

hold on
y = ylim();
y = [y(2),y(2)*.95];
x = bpod.microInjectionStart;
plot([x,x], y,'r')
text(x,y(1), 'inj')
x = bpod.sipperStart;
plot([x,x], y, 'r')
text(x,y(1), 'sip')
x = bpod.tailStart;
plot([x,x], y, 'r')
text(x,y(1), 'tail')
x = bpod.postInjectionStart;
plot([x,x], y, 'r')
text(x,y(1), 'post-inj')

x = bpod.dropDeployed;
y = ones(length(x),1) *y(2); 
scatter(x,y, 'r.')

legend('control','drink high','drink low')
ylabel(yname)
xlabel('time (s)')
xlim([binEdges(1), binEdges(end)])
hold off