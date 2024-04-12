figFolder = 'C:\Users\david\OneDrive - Indiana University\localAlcohol\Figures\2_Sipper\matlabExports\';
% figFolder = 'C:\Users\dis006\OneDrive - Indiana University\localAlcohol\Figures\2_Sipper\matlabExports\';
load("goodClusters.mat")
load("group.mat")

controlColor = [64, 4, 86]/255;
drinkColor = [3, 104, 67]/255;
injectColor = [217,95,2]/255;

isControl = contains(goodClusters.group,'control');
isDrink = contains(goodClusters.group,'drink');
isInject = contains(goodClusters.group,'inject');

%% Fluid consumed

group.mlPerkg = 1000 * group.fluidConsumed ./ group.mouseWeight;

controlConsumed = group.mlPerkg(contains(group.group,'control'));
injectedConsumed = group.mlPerkg(contains(group.group, 'injected'));
drinkConsumed = group.mlPerkg(contains(group.group, 'drink'));


figure(123);clf; hold on
groupNames = ["Control", "Drink", "Inject"];
x = categorical(groupNames);
x = reordercats(x,groupNames);
y = [mean(controlConsumed), mean(drinkConsumed), mean(injectedConsumed)];
err = [std(controlConsumed)/sqrt(length(controlConsumed)), std(drinkConsumed)/sqrt(length(drinkConsumed)), std(injectedConsumed)/sqrt(length(injectedConsumed))];

b = bar(x,y);
b.FaceColor = 'flat';
b.CData(1,:) = controlColor;
b.CData(2,:) = drinkColor;
b.CData(3,:) = injectColor;
errorbar(x,y,err,'k', 'LineStyle','none')
% 
scatter(ones(length(controlConsumed),1)*.75 + .5*rand(length(controlConsumed),1) , controlConsumed, 'k.')
scatter(ones(length(drinkConsumed),1)*1.75 + .5*rand(length(drinkConsumed),1) , drinkConsumed, 'k.')
scatter(ones(length(injectedConsumed),1)*2.75 + .5*rand(length(injectedConsumed),1) , injectedConsumed, 'k.')

ylabel('Fluid consumed (ml/kg)')

f.Units = "inches";
f.Position = [2,2,2,2];
exportgraphics(gcf,[figFolder,'volumeConsumed.pdf'],"BackgroundColor","none","ContentType","vector")

%% brain ethanol 
drinkEthanol = group.brainAlcohol(group.brainAlcohol>0); %hacky because they are the only ones with detectable ethanol

figure(123);clf; hold on
groupNames = ["Control", "Drink", "Inject"];
x = categorical(groupNames);
x = reordercats(x,groupNames);
y = [0, mean(drinkEthanol), 0];
err = [0, std(drinkEthanol)/sqrt(length(drinkEthanol)), 0];

b = bar(x,y);
b.FaceColor = 'flat';
b.CData(1,:) = controlColor;
b.CData(2,:) = drinkColor;
b.CData(3,:) = injectColor;
errorbar(x,y,err,'k', 'LineStyle','none')
% 

scatter(ones(length(drinkEthanol),1)*1.75 + .5*rand(length(drinkEthanol),1) , drinkEthanol, 'k.')


ylabel('Brain [EtOH] (mg/dL)')

f.Units = "inches";
f.Position = [2,2,2,2];
exportgraphics(gcf,[figFolder,'brainEthanol.pdf'],"BackgroundColor","none","ContentType","vector")
%% spikes around drinking 
binWidth=10;
%binEdges = (-60*10):binWidth:(60*12);
target = 'sipperStart';
binEdges = -200:binWidth:1200;
% binEdges = 0:binWidth:60*12;
cLabel = 'zscore';
cRange = [-.4,6];

[spkCounts,  bpod]  = binAroundTarget(goodClusters, target, binEdges);%, 'smooth');
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
yline(0)

plot([0,900], [1,1], 'Color',[.5,.5,.5], 'LineWidth',2)
text(450,1,"Sipper active","HorizontalAlignment","center","VerticalAlignment","bottom")

xlim([x_time(1), x_time(end)])
ylim([-.8,1.1])

xlabel('Time (s)')
ylabel('Z-score')

leg = legend('Control','Drink','Inject','Location','northeast');
legend('boxoff')
leg.ItemTokenSize = [5,4];

f.Units = "inches";
f.Position = [2,2,4.1,2];
exportgraphics(gcf,[figFolder,'grandMean.pdf'], "ContentType","vector","BackgroundColor","none")

%% Spikes around sipper valve
binWidth=0.1;
binEdges = -1:binWidth:14;
spkCounts  = binAroundValve(goodClusters, binEdges);%,'smooth');
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
% title('Average spiking after valve open')
hold on
x_time = binEdges(1:end-1) + (binEdges(2)-binEdges(1))/2;
addShadedLine(x_time, control, {'Color', controlColor})
addShadedLine(x_time, drink, {'Color', drinkColor})
addShadedLine(x_time, inject, {'Color', injectColor})
yline(0)

ylim([-1.1,2])
plot([0,0],ylim,'--','Color',[.5,.5,.5])
plot([10.1,10.1],ylim,'--','Color',[.5,.5,.5])
text(0.1,1.7,["Fluid" "deployed"])
text(10.2,1.7,["Fluid" "removed"])

ylabel('Z-score')
xlabel('Time (s)')
xlim([x_time(1), x_time(end)])


leg = legend('Control', 'Drink', 'Inject',Location='north');
legend('boxoff')
leg.ItemTokenSize = [5,4];

f = gcf;
f.Units = "inches";
f.Position = [2,2,3.3,1.75];
exportgraphics(gcf,[figFolder,'valveMean.pdf'],"ContentType","vector","BackgroundColor","none")

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

%% Spikes - grouped by valve #
figure(8); clf
t = tiledlayout(3,1,'TileSpacing','Compact','Padding','Compact');

nexttile
hold on
y = spkZ_clumped(isControl,:,1);
addShadedLine(x_time, y, {'Color',controlColor})
y = spkZ_clumped(isDrink,:,1);
addShadedLine(x_time, y, {'Color', drinkColor})
y = spkZ_clumped(isInject,:,1);
addShadedLine(x_time, y, {'Color', injectColor})

yline(0)
ylabel('Z-score')
plot([0,0],ylim,'--','Color',[.5,.5,.5])
plot([10.1,10.1],ylim,'--','Color',[.5,.5,.5])
xlim([x_time(1), x_time(end)])
% 
% title('First 12 valve openings')

nexttile
hold on
y = spkZ_clumped(isControl,:,3);
addShadedLine(x_time, y, {'Color',controlColor})
y = spkZ_clumped(isDrink,:,3);
addShadedLine(x_time, y, {'Color',drinkColor})
y = spkZ_clumped(isInject,:,3);
addShadedLine(x_time, y, {'Color',injectColor})

yline(0)
plot([0,0],ylim,'--','Color',[.5,.5,.5])
plot([10.1,10.1],ylim,'--','Color',[.5,.5,.5])
xlim([x_time(1), x_time(end)])
ylabel('zscore')

nexttile
hold on
y = spkZ_clumped(isControl,:,5);
addShadedLine(x_time, y,  {'Color',controlColor})
y = spkZ_clumped(isDrink,:,5);
addShadedLine(x_time, y, {'Color',drinkColor})
y = spkZ_clumped(isInject,:,5);
addShadedLine(x_time, y, {'Color',injectColor})

yline(0)
plot([0,0],ylim,'--','Color',[.5,.5,.5])
plot([10.1,10.1],ylim,'--','Color',[.5,.5,.5])
xlim([x_time(1), x_time(end)])
ylabel('zscore')


%legend('Control', 'Drink', 'Inject')

xlabel('Time (s)')

% leg = legend('Control', 'Drink', 'Inject');
% legend('boxoff')
% leg.ItemTokenSize = [5,4];

% f = gcf;
% f.Units = "inches";
% f.Position = [2,2,3.3,2];
% exportgraphics(gcf,[figFolder,'clumpedMean.pdf'],"ContentType","vector","BackgroundColor","none")

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
xlim([.50,5.5])

f = gcf;
f.Units = "inches";
f.Position = [2,2,1.1,1.7];
exportgraphics(gcf,[figFolder,'pca_scree.pdf'],"ContentType","vector","BackgroundColor","none")

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
yline(0)

ylim([-15,30])
plot([0,0],ylim,'--','Color',[.5,.5,.5])
plot([10.1,10.1],ylim,'--','Color',[.5,.5,.5])

leg = legend({'PC1';'PC2';'PC3'},'Location','best');
legend('boxoff')
leg.ItemTokenSize = [5,4];
xlim([binEdges(1), binEdges(end)])


f = gcf;
f.Units = "inches";
f.Position = [2,2,3.3,1.75];
exportgraphics(gcf,[figFolder,'pca_score.pdf'], "ContentType","vector", "BackgroundColor","none")



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

xline(0)
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
f.Position = [2,2,2,1.7];
exportgraphics(gcf,[figFolder,'pca_PC1Mountain.pdf'])

%%
function vals = val2Ind(x, vals)
    for i=1:length(vals)
        [~,minInd]=min(abs(x-vals(i)));
        vals(i)=minInd;
    end
end

