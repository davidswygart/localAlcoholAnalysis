figFolder = 'C:\Users\david\OneDrive - Indiana University\localAlcohol\Figures\2_Sipper\matlabExports\';
figFolder = 'C:\Users\dis006\OneDrive - Indiana University\localAlcohol\Figures\2_Sipper\matlabExports\';
load("goodClusters.mat")

controlColor = [64, 4, 86]/255;
drinkColor = [3, 104, 67]/255;
injectColor = [217,95,2]/255;

isControl = contains(goodClusters.group,'control');
isDrink = contains(goodClusters.group,'drink');
isInject = contains(goodClusters.group,'inject');
%% spikes around drinking 
binWidth=10;
%binEdges = (-60*10):binWidth:(60*12);
target = 'sipperStart';
binEdges = -4*60:binWidth:60*20;
% binEdges = 0:binWidth:60*12;
cLabel = 'zscore';
cRange = [-.4,6];

[spkCounts,  bpod]  = binAroundTarget(goodClusters, target, binEdges, 'smooth');
spkZ = zscore(spkCounts,0, 2);

figure(1); clf
plotSpikeHeatmapWithBpod(spkZ, binEdges, bpod, cLabel,cRange)
xlabel('Time (10s bins)')

%% spikes around drinking (average)
f = figure(123); clf
hold on
x_time = binEdges(1:end-1) + (binEdges(2)-binEdges(1))/2;
addShadedLine(x_time,spkZ(isControl,:),{'Color', controlColor});
addShadedLine(x_time,spkZ(isDrink,:),{'Color', drinkColor});
addShadedLine(x_time,spkZ(isInject,:),{'Color', injectColor});
plot(xlim(), [0,0], '--k')

scatter(bpod.dropDeployed,ones(length(bpod.dropDeployed)), '.', 'MarkerEdgeColor',[.5,.5,.5])

xlim([x_time(1), x_time(end)])
ylim([-.8,1])

xlabel('Time (s)')
ylabel('Z-score')

leg = legend('Control','Drink','Inject','Location','eastoutside');
legend('boxoff')
leg.ItemTokenSize = [5,4];

f.Units = "inches";
f.Position = [2,2,4,1.75];
exportgraphics(gcf,[figFolder,'grandMean.pdf'])

%% Spikes around sipper valve
binWidth=0.1;
binEdges = -1:binWidth:14;
spkCounts  = binAroundValve(goodClusters, binEdges,'smooth');
spkAvg = mean(spkCounts ,3);
spkZ = zscore(spkAvg,0, 2);


control = spkZ(isControl,:);
drink = spkZ(isDrink,:);
inject = spkZ(isInject,:);

%% Spikes around sipper valve - heatmap
figure(1); clf
t = tiledlayout(3,1,'TileSpacing','Compact','Padding','Compact');

cLabel = 'zscore';
cRange = [-.5,8];

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

%% Spikes around sipper valve - mean
figure(123); clf
title('Average spiking after valve open')
hold on
x_time = binEdges(1:end-1) + (binEdges(2)-binEdges(1))/2;
addShadedLine(x_time, control, {'Color', controlColor})
addShadedLine(x_time, drink, {'Color', drinkColor})
addShadedLine(x_time, inject, {'Color', injectColor})
plot(xlim(), [0,0], '--k')


ylabel('Z-score')
xlabel('Time (s)')


leg = legend('Control', 'Drink', 'Inject');
legend('boxoff')
leg.ItemTokenSize = [5,4];

f = gcf;
f.Units = "inches";
f.Position = [2,2,3,1.75];
exportgraphics(gcf,[figFolder,'valveMean.pdf'])

%% Group by N valves 
valveBinning = 12;
clumped = nan(size(spkCounts,1),size(spkCounts,2),size(spkCounts,3)/valveBinning);
for i=1:size(spkCounts,3)/valveBinning
    ind1 = (i-1)*valveBinning+1;
    ind2 = ind1 + valveBinning - 1;
    clumped(:,:,i) = mean(spkCounts(:,:,ind1:ind2), 3);
end

spkZ_clumped = zscore(clumped,0,2);

% Plot traces for each valve opening with clumped valves
controlValves = spkZ_clumped(isControl, :, :);
controlValves = squeeze(mean(controlValves,1))';

drinkValves = spkZ_clumped(isDrink, :, :);
drinkValves = squeeze(mean(drinkValves,1))';

injectValves = spkZ_clumped(isInject, :, :);
injectValves = squeeze(mean(injectValves,1))';

figure(6); clf
t = tiledlayout(3,1,'TileSpacing','Compact','Padding','Compact');

cLabel = 'Z-score';
cRange = [-1,1];

nexttile
plotSpikeHeatmapWithBpod(controlValves, binEdges, [], cLabel,cRange)
x = [(0-binEdges(1))/binWidth, (10.1-binEdges(1))/binWidth];
scatter(x, [0,0], 'r*')
yL = ylim();
ylim([0,yL(2)])
set(gca,'xticklabel',{[]})
title('Control')
ylabel('Valve number')

nexttile
plotSpikeHeatmapWithBpod(drinkValves, binEdges, [], cLabel,cRange)
x = [(0-binEdges(1))/binWidth, (10.1-binEdges(1))/binWidth];
scatter(x, [0,0], 'r*')
yL = ylim();
ylim([0,yL(2)])
set(gca,'xticklabel',{[]})
title('Drink')
ylabel('Valve number')

nexttile
plotSpikeHeatmapWithBpod(injectValves, binEdges, [], cLabel,cRange)
x = [(0-binEdges(1))/binWidth, (10.1-binEdges(1))/binWidth];
scatter(x, [0,0], 'r*')
yL = ylim();
ylim([0,yL(2)])
set(gca,'xticklabel',{[]})
title('Inject')
ylabel('Valve number')
xlabel('Time (0.1s bins)')

%% Pull out peak Inh
%inhWindow = val2Ind(x_time, [1.6,3.2]); %time where zscore is <-0.5 for any trace in Fig 4
inhWindow = val2Ind(x_time, [1.25,3.55]); %time below 0
inh = squeeze(mean(spkZ_clumped(:,inhWindow(1):inhWindow(2),:), 2));


figure(7);clf

hold on
title('"Inh" peak (1.6-3.2s after valve)')
addShadedLine([],inh(isControl,:), {'Color', controlColor})
addShadedLine([],inh(isDrink,:),{'Color', drinkColor})
addShadedLine([],inh(isInject,:),{'Color', injectColor})
plot(xlim(), [0,0], '--k')
legend('control', 'drink', 'inject')
xlabel('Time (valve group #)')
ylabel('zscore')
hold off

%% Spikes around sipper valve - mean
figure(123); clf
% t = tiledlayout(3,1,'TileSpacing','Compact','Padding','Compact');
% 
% nexttile
hold on
y = spkZ_clumped(isControl,:,1);
addShadedLine(x_time, y, {'Color',controlColor})
y = spkZ_clumped(isDrink,:,1);
addShadedLine(x_time, y, {'Color', drinkColor})
y = spkZ_clumped(isInject,:,1);
addShadedLine(x_time, y, {'Color', injectColor})
plot(xlim(), [0,0], '--k')
ylabel('Z-score')

title('First 12 valve openings')

% nexttile
% hold on
% y = spkZ_clumped(isControl,:,3);
% addShadedLine(x_time, y, 'b','Control')
% y = spkZ_clumped(isDrink,:,3);
% addShadedLine(x_time, y, 'r', 'Drink')
% y = spkZ_clumped(isInject,:,3);
% addShadedLine(x_time, y, 'g', 'Inject')
% plot(xlim(), [0,0], '--k')
% ylabel('zscore')
% 
% nexttile
% hold on
% y = spkZ_clumped(isControl,:,5);
% addShadedLine(x_time, y, 'b','Control')
% y = spkZ_clumped(isDrink,:,5);
% addShadedLine(x_time, y, 'r', 'Drink')
% y = spkZ_clumped(isInject,:,5);
% addShadedLine(x_time, y, 'g', 'Inject')
% plot(xlim(), [0,0], '--k')
% ylabel('zscore')


%legend('Control', 'Drink', 'Inject')

xlabel('Time (s)')

% leg = legend('Control', 'Drink', 'Inject');
% legend('boxoff')
% leg.ItemTokenSize = [5,4];

f = gcf;
f.Units = "inches";
f.Position = [2,2,3,2];
exportgraphics(gcf,[figFolder,'clumpedMean.pdf'])

%% Run PCA
[coeff,score,latent,~,explained] = pca(spkZ'); % <-- Data are (time(bins) x neurons)

%% Plot explained variance or Scree
figure(123);clf
% plot(latent,'.-')
ylabel('Eigenvalue')
bar(explained,'k');  
% bar(cumsum(explained),'k');
ylabel('Explained variance (%)')
xlabel('PC')
xlim([.50,6.5])

f = gcf;
f.Units = "inches";
f.Position = [2,2,2,2];
exportgraphics(gcf,[figFolder,'pca_scree.pdf'])

%% Plot PC pattern (score)
figure(123); clf
hold on

x_time = binEdges(1:end-1) + (binEdges(2)-binEdges(1))/2;
plot(x_time, score(:,1),'LineWidth',2,'Color',[.8,0,.2,1]);
plot(x_time, score(:,2),'LineWidth',2,'Color',[0,0,1,.5]);
plot(x_time, score(:,3),'LineWidth',2,'Color',[0,.8,.2,.5]);
%plot(x_time, score(:,4),'LineWidth',2,'Color',[.5,.5,.5,.5]);

xlabel('Time (s)');
ylabel('Z-score');
xlim([binEdges(1), binEdges(end)])

hold on
plot(xlim,[0,0],'--k')

leg = legend({'PC1';'PC2';'PC3'},'Location','best');
legend('boxoff')
leg.ItemTokenSize = [5,4];
xlim([binEdges(1), binEdges(end)])


f = gcf;
f.Units = "inches";
f.Position = [2,2,3,2];
exportgraphics(gcf,[figFolder,'pca_score.pdf'])

%% Compare PC1 vs PC2, split by type

figure(4);clf
scatter(coeff(isControl,1),coeff(isControl,2), 'b','filled','MarkerFaceAlpha',.5)
hold on
scatter(coeff(isDrink,1),coeff(isDrink,2), 'r','filled','MarkerFaceAlpha',.5)
scatter(coeff(isInject,1),coeff(isInject,2), 'g','filled','MarkerFaceAlpha',.5)


plot([-1,1], [0,0], 'k--')
plot([0,0],[1,-1], 'k--')


xlim([-.17,.17])
ylim([-.17,.17])
axis equal 

xlabel("PC1")
ylabel("PC2"),

legend('Control', 'Drink', 'Inject','Location','eastoutside')

%% Compare PC1 vs PC4, split by type

figure(4);clf
scatter(coeff(isControl,1),coeff(isControl,4), 'b','filled','MarkerFaceAlpha',.5)
hold on
scatter(coeff(isDrink,1),coeff(isDrink,4), 'r','filled','MarkerFaceAlpha',.5)
scatter(coeff(isInject,1),coeff(isInject,4), 'g','filled','MarkerFaceAlpha',.5)

plot([-1,1], [0,0], 'k--')
plot([0,0],[1,-1], 'k--')

xlim([-.17,.17])
ylim([-.17,.17])
axis equal 

xlabel("PC1")
ylabel("PC4"),

legend('Control', 'Drink', 'Inject','Location','eastoutside')

%% PC1 loadings, split by group  (mountain plot)
figure(123); clf
[f,x] = ecdf(coeff(isControl,1));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',1.5,'Color',controlColor)
hold on

[f,x] = ecdf(coeff(isDrink,1));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',1.5,'Color',drinkColor)

[f,x] = ecdf(coeff(isInject,1));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',1.5,'Color', injectColor)

plot([0,0],[0,0.5], '--k')
xlabel('PC1 loading (coeff)')
ylabel('Folded probability')
legend("Control","Drink", "Inject","Location","northwest")

[h,p] = kstest2(coeff(isControl,1),coeff(isDrink,1) );
% text(-.14,.4,['CvD p=',num2str(p*2,2)], 'FontSize',5)

[h,p] = kstest2(coeff(isControl,1),coeff(isInject,1) );
% text(-.14,.35,['CvI p=',num2str(p*2,2)], 'FontSize',5)

% pc1_Thresh = -.05;
% plot([pc1_Thresh,pc1_Thresh], [0,.5], '-.k')
leg = legend("Control","Drink", "Inject","Location","northeast");
legend('boxoff')
leg.ItemTokenSize = [5,4];

f = gcf;
f.Units = "inches";
f.Position = [2,2,2,2];
exportgraphics(gcf,[figFolder,'pca_PC1Mountain.pdf'])

%%
function vals = val2Ind(x, vals)
    for i=1:length(vals)
        [~,minInd]=min(abs(x-vals(i)));
        vals(i)=minInd;
    end
end

