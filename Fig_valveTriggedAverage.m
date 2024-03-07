%load("allClusters.mat")

allClusters = allClusters(allClusters.fr>0.1, :);
allClusters = allClusters(~contains(allClusters.phyLabel,'noise'), :);
allClusters = allClusters(allClusters.presenceRatio>0.9, :);
allClusters = allClusters(allClusters.isiViolations<1, :);

%% choose bin width
binWidth=0.1;
binEdges = -1:binWidth:14;

%%
includedValves = 1:60;
spkCounts  = binAroundValve(allClusters, binEdges,'smooth');

% average and then zscore (control peaks at 0.6 std)
%spkAvg = mean(spkCounts(:,:,includedValves) ,3);
%spkZ = zscore(spkAvg,0, 2);

% zscore and then average (control peaks at 0.1 std)
spkZ = zscore(spkCounts,0,2);
spkZ = mean(spkZ(:,:,includedValves),3);


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
%% Figure 4 - spiking for all valves averaged 
figure(4)
clf
title('Average spiking after valve open')
shadedErrorBar(binEdges(1:end-1), avgControl,semControl, 'lineProps', 'b')
hold on
shadedErrorBar(binEdges(1:end-1), avgDrink,semDrink, 'lineProps', 'r')
shadedErrorBar(binEdges(1:end-1), avgInject,semInject, 'lineProps', 'g')

legend('control','drink','inject')
ylabel(yname)
xlabel('time (s)')
xlim([binEdges(1), binEdges(end)])
hold off

%% save for stats in R
%saveCsvForR_repeatedMeasuresMixedModel({control,drink,inject},{"control","drink","inject"},'averageAroundValve.csv')

%% run multiple ttests
nTimepoints = size(control,2);
pVals = nan(nTimepoints,3);
for i=1:nTimepoints
    [~,c_i]=ttest2(control(:,i), inject(:,i),'Vartype','unequal');
    [~,c_d]=ttest2(control(:,i), drink(:,i),'Vartype','unequal');
    [~,i_d]=ttest2(inject(:,i), drink(:,i),'Vartype','unequal');
    pVals(i,:) = [c_i,c_d,i_d];
end
h = fdr_bh(pVals);

hold on
x = binEdges(1:end-1);
y = ylim();
y = ones(length(x),1) * y(1);
scatter(x(any(h,2)),y(any(h,2)),'*k')
hold off