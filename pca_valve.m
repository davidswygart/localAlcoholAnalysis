%% Spikes around sipper valve
binWidth=0.1;
binEdges = -1:binWidth:14;
spkCounts  = binAroundValve(goodClusters, binEdges,'smooth');
spkAvg = mean(spkCounts ,3);
spkZ = zscore(spkAvg,0, 2);

isControl = contains(goodClusters.group,'control');
control = spkZ(isControl,:);
isDrink = contains(goodClusters.group,'drink');
drink = spkZ(isDrink,:);
isInject = contains(goodClusters.group,'inject');
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
figure(2); clf
title('Average spiking after valve open')
hold on
x_time = binEdges(1:end-1) + (binEdges(2)-binEdges(1))/2;
addShadedLine(x_time, control, 'b','Control')
addShadedLine(x_time, drink, 'r', 'Drink')
addShadedLine(x_time, inject, 'g', 'Inject')
plot(xlim(), [0,0], '--k')

legend('Control', 'Drink', 'Inject')
ylabel('zscore')
xlabel('time (s)')

%% Run PCA
[coeff,score,latent,~,explained] = pca(spkZ'); % <-- Data are (time(bins) x neurons)

%% Plot explained variance or Scree
figure(25);clf
% plot(latent,'.-')
ylabel('Eigenvalue')
bar(explained,'k');  
% bar(cumsum(explained),'k');
ylabel('Explained variance (%)')
xlabel('PC')
xlim([.50,10.5])

%% Plot PC pattern (score)
figure(3); clf
hold on

x_time = binEdges(1:end-1) + (binEdges(2)-binEdges(1))/2;
plot(x_time, score(:,1),'LineWidth',2,'Color',[.8,0,.2,.5]);
plot(x_time, score(:,2),'LineWidth',2,'Color',[0,0,1,.5]);
plot(x_time, score(:,3),'LineWidth',2,'Color',[0,.8,.2,.5]);
plot(x_time, score(:,4),'LineWidth',2,'Color',[.5,.5,.5,.5]);

xlabel('Time (s)'); ylabel('Z-score');

hold on

legend({'PC1';'PC2';'PC3';'PC4'},'Location','eastoutside')
xlim([binEdges(1), binEdges(end)])

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
figure(5); clf
[f,x] = ecdf(coeff(isControl,1));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',2,'Color','b')
hold on

[f,x] = ecdf(coeff(isDrink,1));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',2,'Color','r')

[f,x] = ecdf(coeff(isInject,1));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',2,'Color','g')

plot([0,0],[0,0.5], '-k')
xlabel('PC1 loading (coeff)')
ylabel('Folded probability')
legend("Control","Drink", "Inject","Location","northwest")

[h,p] = kstest2(coeff(isControl,1),coeff(isDrink,1) );
text(-.14,.35,['CvD p=',num2str(p*2,2)])

[h,p] = kstest2(coeff(isControl,1),coeff(isInject,1) );
text(-.14,.32,['CvI p=',num2str(p*2,2)])

pc1_Thresh = -.05;
plot([pc1_Thresh,pc1_Thresh], [0,.5], '-.k')
legend("Control","Drink", "Inject","Location","northwest")

%% PC2 loadings, split by group  (mountain plot)
figure(5); clf
[f,x] = ecdf(coeff(isControl,2));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',2,'Color','b')
hold on

[f,x] = ecdf(coeff(isDrink,2));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',2,'Color','r')

[f,x] = ecdf(coeff(isInject,2));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',2,'Color','g')

plot([0,0],[0,0.5], '--k')
xlabel('PC2 loading (coeff)')
ylabel('Folded probability')
legend("Control","Drink", "Inject","Location","northwest")

[h,p] = kstest2(coeff(isControl,2),coeff(isDrink,2) );
text(-.09,.35,['CvD p=',num2str(p*2,2)])

[h,p] = kstest2(coeff(isControl,2),coeff(isInject,2) );
text(-.09,.32,['CvI p=',num2str(p*2,2)])

%% PC3 loadings, split by group  (mountain plot)
figure(5); clf
[f,x] = ecdf(coeff(isControl,3));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',2,'Color','b')
hold on

[f,x] = ecdf(coeff(isDrink,3));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',2,'Color','r')

[f,x] = ecdf(coeff(isInject,3));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',2,'Color','g')

plot([0,0],[0,0.5], '--k')
xlabel('PC3 loading (coeff)')
ylabel('Folded probability')
legend("Control","Drink", "Inject","Location","northwest")

[h,p] = kstest2(coeff(isControl,3),coeff(isDrink,3) );
text(-.13,.35,['CvD p=',num2str(p*2,2)])

[h,p] = kstest2(coeff(isControl,3),coeff(isInject,3) );
text(-.13,.32,['CvI p=',num2str(p*2,2)])

%% PC4 loadings, split by group  (mountain plot)
figure(5); clf
[f,x] = ecdf(coeff(isControl,4));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',2,'Color','b')
hold on

[f,x] = ecdf(coeff(isDrink,4));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',2,'Color','r')

[f,x] = ecdf(coeff(isInject,4));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',2,'Color','g')

plot([0,0],[0,0.5], '--k')
xlabel('PC4 loading (coeff)')
ylabel('Folded probability')
legend("Control","Drink", "Inject","Location","northwest")

[h,p] = kstest2(coeff(isControl,4),coeff(isDrink,4) );
text(-.13,.35,['CvD p=',num2str(p*2,2)])

[h,p] = kstest2(coeff(isControl,4),coeff(isInject,4) );
text(-.13,.32,['CvI p=',num2str(p*2,2)])

%% Plot traces for each valve opening, only for heavy negative PC1 loaders
spkZ_valves = zscore(spkCounts,0,2);

controlValves = spkZ_valves(isControl & (coeff(:,1) < pc1_Thresh), :, :);
controlValves = squeeze(mean(controlValves,1))';

drinkValves = spkZ_valves(isDrink & (coeff(:,1) < pc1_Thresh), :, :);
drinkValves = squeeze(mean(drinkValves,1))';

injectValves = spkZ_valves(isInject & (coeff(:,1) < pc1_Thresh), :, :);
injectValves = squeeze(mean(injectValves,1))';

figure(6); clf
t = tiledlayout(3,1,'TileSpacing','Compact','Padding','Compact');

cLabel = 'zscore';
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

%% Group by N valves 
valveBinning = 10;
clumped = nan(size(spkCounts,1),size(spkCounts,2),size(spkCounts,3)/valveBinning);
for i=1:size(spkCounts,3)/valveBinning
    ind1 = (i-1)*valveBinning+1;
    ind2 = ind1 + valveBinning - 1;
    clumped(:,:,i) = mean(spkCounts(:,:,ind1:ind2), 3);
end

spkZ_clumped = zscore(clumped,0,2);

%% Plot traces for each valve opening, only heavy negative PC1 loaders, with clumped valves
controlValves = spkZ_clumped(isControl & (coeff(:,1) < pc1_Thresh), :, :);
controlValves = squeeze(mean(controlValves,1))';

drinkValves = spkZ_clumped(isDrink & (coeff(:,1) < pc1_Thresh), :, :);
drinkValves = squeeze(mean(drinkValves,1))';

injectValves = spkZ_clumped(isInject & (coeff(:,1) < pc1_Thresh), :, :);
injectValves = squeeze(mean(injectValves,1))';

figure(6); clf
t = tiledlayout(3,1,'TileSpacing','Compact','Padding','Compact');

cLabel = 'zscore';
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

%% Pull out peak Inh and Exc Values for heavy negative loaders
%inhWindow = val2Ind(x_time, [1.6,3.2]); %time where zscore is <-0.5 for any trace in Fig 4
inhWindow = val2Ind(x_time, [1.25,3.55]); %time below 0
inh = squeeze(mean(spkZ_clumped(:,inhWindow(1):inhWindow(2),:), 2));


figure(7);clf

hold on
title('"Inh" peak (1.6-3.2s after valve), heavy negative loaders only')
addShadedLine([],inh(isControl & (coeff(:,1) < pc1_Thresh),:),'b','Control')
addShadedLine([],inh(isDrink & (coeff(:,1) < pc1_Thresh),:),'r','Drink')
addShadedLine([],inh(isInject & (coeff(:,1) < pc1_Thresh),:),'g','Inject')
plot(xlim(), [0,0], '--k')
legend('control', 'drink', 'inject')
xlabel('Time (valve group #)')
ylabel('zscore')
hold off



%%
function vals = val2Ind(x, vals)
    for i=1:length(vals)
        [~,minInd]=min(abs(x-vals(i)));
        vals(i)=minInd;
    end
end

