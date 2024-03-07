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

%% Figure 4 - spiking for all valves averaged 
figure(4)
clf
title('Average spiking after valve open')
hold on
x = binEdges(1:end-1);
addShadedLine(x, control, 'b','Control')
addShadedLine(x, drink, 'r', 'Drink')
addShadedLine(x, inject, 'g', 'Inject')

legend()
ylabel(yname)
xlabel('time (s)')

%plot statistically significant points
x = binEdges(1:end-1);
x = x(any(h,2));
y = ylim();
y = ones(length(x),1) * y(1);
scatter(x,y,'*k')
hold off

%% save for stats in R
%saveCsvForR_repeatedMeasuresMixedModel({control,drink,inject},{"control","drink","inject"},'averageAroundValve.csv')

