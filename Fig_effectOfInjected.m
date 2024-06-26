%load("allClusters.mat")
allClusters = allClusters(allClusters.fr>0.1, :);
allClusters = allClusters(~contains(allClusters.phyLabel,'noise'), :);
allClusters = allClusters(allClusters.presenceRatio>0.9, :);
allClusters = allClusters(allClusters.isiViolations<1, :);

%% choose anchor point and bin width
binWidth=10;

% binEdges = -60*10:binWidth:60*25;
% target = 'sipperStart';

binEdges = -200:binWidth:60*12;
target = 'microInjectionStart';
 %% Choose spike rate or z score

[spkCounts,  bpod]  = binAroundTarget(allClusters, target, binEdges, 'smooth');
spkRate = spkCounts / binWidth;

spkZ = zscore(spkCounts,0, 2);
yname = 'zscore';

% isBaseline = binEdges(1:end-1) < 0;
% baseline = spkCounts(:,isBaseline);
% avgBaseline = mean(baseline,2);
% %stdBaseline = std(baseline,0,2);
% spkZ = (spkCounts-avgBaseline)./std(spkCounts,0,2);
% yname = 'zscore (baseline subtracted)';

% spkZ = spkCounts;
% yname = 'spike count';

 %% Choose spike rate or z score
control = spkZ(contains(allClusters.group,'control') | contains(allClusters.group,'drink'), :);
inject = spkZ(contains(allClusters.group,'injected'), :);

avgControl = mean(control,1);
stdControl = std(control,0, 1);
semControl = stdControl / sqrt(size(control,1));

avgInject =  mean(inject,1);
stdInject =  std(inject,0, 1);
semInject = stdInject / sqrt(size(inject,1));

figure(1)
clf
shadedErrorBar(binEdges(1:end-1), avgControl,semControl, 'lineProps', 'b')
hold on
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
scatter(x,y, 'r.')

legend('control and drink','injected')
ylabel(yname)
xlabel('time (s)')
xlim([binEdges(1), binEdges(end)])
hold off
%%
saveCsvForR_repeatedMeasuresMixedModel({control,inject},{"control","inject"}, 'effectOfInject.csv')
%% perform ttests
nTimepoints = size(control,2);
pVals = nan(nTimepoints,1);
for i=1:nTimepoints
    [~,p]=ttest2(control(:,i), inject(:,i),'Vartype','unequal');
    pVals(i) = p;
end
h = fdr_bh(pVals);

hold on
x = binEdges(1:end-1);
y = ylim();
y = ones(length(x),1) * y(1);
scatter(x(h),y(h),'*k')
hold off
