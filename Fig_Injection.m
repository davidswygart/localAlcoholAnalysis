%%
figFolder = 'C:\Users\david\OneDrive - Indiana University\localAlcohol\Figures\1_Injection\matlabExports\';
load("diffusion.mat")
load("goodClusters.mat")
%% spikes around injection 
binWidth=10;
%binEdges = (-60*10):binWidth:(60*12);
target = 'microInjectionStart';
binEdges = -200:binWidth:60*12;
% binEdges = 0:binWidth:60*12;
cLabel = 'zscore';
cRange = [-.4,6];

[spkCounts,  bpod]  = binAroundTarget(goodClusters, target, binEdges, 'smooth');
spkZ = zscore(spkCounts,0, 2);

figure(1); clf
plotSpikeHeatmapWithBpod(spkZ, binEdges, bpod, cLabel,cRange)
xlabel('Time (10s bins)')

isControl = contains(goodClusters.group,'control') | contains(goodClusters.group,'drink');
isInject = contains(goodClusters.group,'inject');

%% spikes around injection (average)
f = figure(123); clf
x_time = binEdges(1:end-1) + (binEdges(2)-binEdges(1))/2;
addShadedLine(x_time,spkZ(isControl,:),'b','Control');
hold on
addShadedLine(x_time,spkZ(isInject,:),'r','Inject');
plot(xlim(), [0,0], '--k')

%plotBpod(bpod)
xlim([binEdges(1), binEdges(end)])

xlabel('Time (s)')
ylabel('zscore')

leg = legend('Control','Inject','Location','best');
legend('boxoff')
leg.ItemTokenSize = [5,4];

f.Units = "inches";
f.Position = [2,2,3,2];
exportgraphics(gcf,[figFolder,'grandMean.pdf'])

%% Run PCA
[coeff,score,latent,~,explained] = pca(spkZ'); % <-- Data are (time(bins) x neurons)

%% Plot explained variance or Scree
f = figure(123);clf
% plot(latent,'.-')
ylabel('Eigenvalue')
bar(explained,'k');  
% bar(cumsum(explained),'k');
ylabel('Explained variance (%)')
xlabel('PC')
xlim([.50,6.5])

f.Units = "inches";
f.Position = [2,2,2,2];
exportgraphics(gcf,[figFolder,'Scree.pdf'])

%% Plot PC pattern (score)
f = figure(123);clf
hold on

x_time = binEdges(1:end-1) + (binEdges(2)-binEdges(1))/2;
plot(x_time, score(:,1),'LineWidth',1.5,'Color',[.8,0,.2,.3]);
plot(x_time, score(:,2),'LineWidth',2,'Color',[0,0,1,1]);
%plot(x_time, score(:,3),'LineWidth',1.5,'Color',[0,.8,.2,.3]);

xlabel('Time (s)'); ylabel('Z-score');

hold on
plot(xlim,[0,0], '--k')
%plotBpod(bpod)

leg = legend({'PC1';'PC2'},'Location','best');
legend('boxoff')
leg.ItemTokenSize = [5,4];
xlim([binEdges(1), binEdges(end)])

f.Units = "inches";
f.Position = [2,2,3,2];
exportgraphics(gcf,[figFolder,'pca_score.pdf'])

%% PC2 loadings, split by group  (mountain plot)
figure(123);clf
[f,x] = ecdf(coeff(isControl,2));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',2,'Color','b')
hold on
[f,x] = ecdf(coeff(isInject,2));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',2,'Color','r')

plot([0,0],[0,0.5], '--k')
xlabel('PC2 loading (coeff)')
ylabel('Folded probability')
leg = legend("Control", "Inject","Location","northeast");
legend('boxoff')
leg.ItemTokenSize = [5,4];

[h,p] = kstest2(coeff(isControl,2),coeff(isInject,2) );
text(-.1,.4,['p=',num2str(p,2)], 'FontSize',5)
f = gcf;
f.Units = "inches";
f.Position = [2,2,2,2];
exportgraphics(gcf,[figFolder,'pca_loading.pdf'])


%% plot concentration at different distances
figure(123); clf
hold on
thisDist = 39;%385 um
plot(diffusion.time, diffusion.conc(thisDist,:))
thisDist = 56;%560 um
plot(diffusion.time, diffusion.conc(thisDist,:))
thisDist = 78;%775 um
plot(diffusion.time, diffusion.conc(thisDist,:))
plot([120,120],ylim,'r--')
xlabel('Time (s)')
ylabel('Concentration (mg/dL)')

leg = legend('380 um', '560 um', '775 um');
legend('boxoff')
leg.ItemTokenSize = [5,4];
f = gcf;
f.Units = "inches";
f.Position = [2,2,2,2];
exportgraphics(gcf,[figFolder,'concentration_atDistances.pdf'])

%% Calculate predicted alcohol time per cluster
conc = calcConcentrationForCluster(goodClusters, diffusion, binEdges);

figure(3); clf

plotSpikeHeatmapWithBpod(conc, binEdges, bpod, 'Conentration (mg/dL)' ,[0,300])
xlabel('Time (10s bins)')
%% scatter plot: alcohol vs. spiking
figure(5); clf
x = conc(isInject,:);
y = spkZ(isInject,:);
scatter(x(:),y(:),5,'filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',0)

xlabel('Concentration (mg/dL)')
ylabel('Spiking (Z-score)')

[p,S] = polyfit(x(:),y(:),2);
xp = 0:2:160;
[yp, delta] = polyval(p, xp, S);
hold on
shadedErrorBar(xp,yp,delta)
xlim([0,165])


%% Correlation of firing to concentration (lag = 0)
 [rho,pval] = corr(conc',spkZ');

 identity = logical(eye(length(rho)));
 pval = pval(identity);
 rho = rho(identity);
%% Mountain plot correlation lag 0
figure(5); clf
[f,x] = ecdf(rho(isControl));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',2,'Color','b')
hold on

[f,x] = ecdf(rho(isInject));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',2,'Color','r')

xlabel('Correlation (firing to ethanol concentration)')
ylabel('Folded probability')

plot([0,0],[0,.5], 'k--')
legend("Control","Inject")

%% Corelation to firing (stacked bar) for paper
figure(123); clf
stack = nan(3,2);
stack(1,1) = sum((rho(isInject)<0) & (pval(isInject)<.05));
stack(2,1) = sum(pval(isInject)>.05);
stack(3,1) = sum((rho(isInject)>0) & (pval(isInject)<.05));
stack(:,1) = stack(:,1) / sum(stack(:,1));

stack(1,2) = sum((rho(isControl)<0) & (pval(isControl)<.05));
stack(2,2) = sum(pval(isControl)>.05);
stack(3,2) = sum((rho(isControl)>0) & (pval(isControl)<.05));
stack(:,2) = stack(:,2) / sum(stack(:,2));

groupNames = ["Inject","Control"];
x = categorical(groupNames);
x = reordercats(x, groupNames);
barh(x, stack', 'stacked')
%ylim([0,sum(control_stack)])
xlabel('Proportion of clusters')

c = colororder;
textY = 2.7;
textX = stack(1,2)/2;
text(textX,textY, '-Corr','HorizontalAlignment','center','Color',c(1,:),'FontSize',8)
textX = stack(1,2) + stack(2,2)/2;
text(textX,textY, 'NS','HorizontalAlignment','center','Color',c(2,:),'FontSize',8)
textX = stack(1,2) + stack(2,2) + stack(3,2)/2;
text(textX,textY, '+Corr','HorizontalAlignment','center','Color',c(3,:),'FontSize',8)

%ylim([.25,3])

f = gcf;
%hA = axes(f);

box off
% set(gca, 'XTick', [], 'XTickLabel', []);
% set(hA, 'YTick', [], 'YTickLabel', []);
% set(get(gca, 'XAxis'), 'Visible', 'off');
% set(get(gca, 'YAxis'), 'Visible', 'off');

f.Units = "inches";
f.Position = [2,2,2,1.5];
exportgraphics(gcf,[figFolder,'concentration_stackedBar.pdf'])


%% Scatter plot correlation vs distance
figure(7); clf

subplot(2,1,1)
x = goodClusters.distFromInj(isControl);
y = rho(isControl);
scatter(x,y, 'filled', 'MarkerEdgeAlpha',0, 'MarkerFaceAlpha',.8)
xlim([200,1000])
title('Control')
ylabel('Correlation: Spiking vs concentration')

subplot(2,1,2)
x = goodClusters.distFromInj(isInject);
y = rho(isInject);
scatter(x,y, 'filled', 'MarkerEdgeAlpha',0, 'MarkerFaceAlpha',.8)
xlim([200,1000])
title('Inject')
ylabel('Correlation: spiking vs concentration')
xlabel('Distance from injection (um)')

%% Cross correlation
[r,lag] = xcorr(conc(1,:), spkZ(1,:));
r = nan(size(conc,1), length(r));
r_shuff = r;

for i=1:size(conc,1)
    r(i,:) = xcorr(conc(i,:), spkZ(i,:), 'normalized');
    spk_shuff = spkZ(i,randperm(size(spkZ,2)));
    r_shuff(i,:) = xcorr(conc(i,:), spk_shuff, 'normalized');

end

%% Plot mean cross correlation (all)
figure(8); clf
x = lag*binWidth;
addShadedLine(x, r, 'k','All')
hold on
addShadedLine(x, r_shuff, 'k--', 'Shuffled')
% ylim([0,1])
plot([0,0],ylim,'k--')
plot(xlim,[0,0], 'k--')
xlabel('Lag (s)')
ylabel('Correlation: spiking vs concentration')

%% Plot mean cross correlation (split by type)
figure(8); clf
x = lag*binWidth;
addShadedLine(x, r(isControl,:), 'b','Control')
hold on
addShadedLine(x, r(isInject,:), 'r','Inject')

addShadedLine(x, r_shuff(isControl,:), 'b--','Control')
addShadedLine(x, r_shuff(isInject,:), 'r--','Inject')

%ylim([0,1])
plot([0,0],ylim,'k--')
plot(xlim,[0,0], 'k--')
xlabel('Lag (s)')
ylabel('Correlation: spiking vs concentration')
legend("Control", "Inject")

%% Mountain plot of p-Value for specific lag
lagPoint = -60;
%lagPoint = -360;
[~,lagInd] = min(abs(lag-lagPoint/binWidth));

figure(9); clf
[f,x] = ecdf(r(isControl,lagInd));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',2,'Color','b')
hold on

[f,x] = ecdf(r(isInject,lagInd));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',2,'Color','r')

xlabel('Correlation (firing to ethanol concentration)')
ylabel('Folded probability')

title(['Spike lag ', num2str(lagPoint), 's'])

plot([0,0],[0,.5], 'k--')
legend("Control","Inject")

%% Scatter plot correlation vs distance for lag point
figure(10); clf

subplot(2,1,1)
x = goodClusters.distFromInj(isControl);
y = r(isControl, lagInd);
scatter(x,y, 'filled', 'MarkerEdgeAlpha',0, 'MarkerFaceAlpha',.8)
xlim([600,1500])
title('Control')
ylabel('Correlation: Spiking vs concentration')

subplot(2,1,2)
x = goodClusters.distFromInj(isInject);
y = r(isInject, lagInd);
scatter(x,y, 'filled', 'MarkerEdgeAlpha',0, 'MarkerFaceAlpha',.8)
xlim([600,1500])
title('Inject')
ylabel('Correlation: spiking vs concentration')
xlabel('Distance from injection (um)')
%% calc concentration function
function conc = calcConcentrationForCluster(clusters, diffusion,binEdges)
    t = binEdges(1:end-1) + (binEdges(2)-binEdges(1))/2;
    
    interpConc = nan(length(diffusion.dist), length(t));
    for i = 1:length(diffusion.dist)
        interpConc(i,:) = interp1(diffusion.time,diffusion.conc(i,:),t,'nearest','extrap');
    end
    
    conc = nan(size(clusters,1), length(t));
    for i = 1:length(clusters.distFromInj)
        clusterDist = clusters.distFromInj(i);
        [~,distInd] = min(abs(diffusion.dist-clusterDist));
        conc(i,:) = interpConc(distInd,:);
    end
end