allClusters = allClusters(allClusters.fr>0.1, :);
allClusters = allClusters(~contains(allClusters.phyLabel,'noise'), :);
allClusters = allClusters(allClusters.presenceRatio>0.9, :);
allClusters = allClusters(allClusters.isiViolations<1, :);

%% choose anchor point and bin width
binWidth=10;

binEdges = -60*17:binWidth:60*25;
target = 'sipperStart';

% binEdges = -60*5:binWidth:60*35;
% target = 'microInjectionStart';
%%

[spkCounts,  bpod]  = binAroundTarget(allClusters, target, binEdges);
spkRate = spkCounts / (binEdges(2)-binEdges(1));
spkZ = zscore(spkCounts,0, 2);

 %% Choose spike rate or z score
control = spkRate(contains(allClusters.group,'control'), :);
drink = spkRate(contains(allClusters.group,'drink'), :);
inject = spkRate(contains(allClusters.group,'injected'), :);
yname = 'spike rate';

% control = spkZ(contains(allClusters.group,'control'), :);
% drink = spkZ(contains(allClusters.group,'drink'), :);
% inject = spkZ(contains(allClusters.group,'injected'), :);
% yname = 'zscore';


%%
avgControl = mean(control,1);
stdControl = std(control,0, 1);
semControl = stdControl / sqrt(size(control,1));

avgDrink =  mean(drink,1);
stdDrink =  std(drink,0, 1);
semDrink = stdDrink / sqrt(size(drink,1));

avgInject =  mean(inject,1);
stdInject =  std(inject,0, 1);
semInject = stdInject / sqrt(size(inject,1));


shadedErrorBar(binEdges(1:end-1), avgControl,semControl, 'lineProps', 'b')
hold on
shadedErrorBar(binEdges(1:end-1), avgDrink,semDrink, 'lineProps', 'g')
shadedErrorBar(binEdges(1:end-1), avgInject,semInject, 'lineProps', 'r')

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
scatter(x,y, 'r*')

legend('control','drink','injected')
ylabel(yname)
xlabel('time (s)')
xlim([binEdges(1), binEdges(end)])
hold off