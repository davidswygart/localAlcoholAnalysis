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
figure(2); clf
x_time = binEdges(1:end-1) + (binEdges(2)-binEdges(1))/2;
addShadedLine(x_time,spkZ(isControl,:),'b','Control')
hold on
addShadedLine(x_time,spkZ(isInject,:),'r','Inject')
plot(xlim(), [0,0], '--k')

plotBpod(bpod)
xlim([binEdges(1), binEdges(end)])

xlabel('time (s)')
ylabel('zscore')

legend('Control','Inject')
%% Calculate predicted alcohol time per cluster
conc = calcConcentrationForCluster(goodClusters, diffusion, binEdges);

figure(3); clf

plotSpikeHeatmapWithBpod(conc, binEdges, bpod, 'Conentration (mg/dL)' ,[0,150])
xlabel('Time (10s bins)')
%% What was the peak and when was it reached
figure(4);clf
[peak, peakInd] = max(conc(isInject,:),[],2);
scatterhist(x_time(peakInd),peak, 'Marker','.')
xlabel('Time of peak concentration (s)')
ylabel('Peak concentration (mg/dL)')

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

%% scatter plot: alcohol vs. spiking (only peak concentration)
[peak, peakInd] = max(conc,[],2);
spkAtPeak = spkZ(:,peakInd);


figure(5); clf
x = peak(isInject);
y = spkAtPeak(isInject);
scatter(x(:),y(:),5,'filled','MarkerFaceAlpha',.7,'MarkerEdgeAlpha',0)

xlabel(' Peak concentration (mg/dL)')
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

%% Corelation to firing (stacked bar)
figure(6); clf
stack = nan(3,2);
stack(1,1) = sum((rho(isControl)<0) & (pval(isControl)<.05));
stack(2,1) = sum(pval(isControl)>.05);
stack(3,1) = sum((rho(isControl)>0) & (pval(isControl)<.05));
stack(:,1) = stack(:,1) / sum(stack(:,1));

stack(1,2) = sum((rho(isInject)<0) & (pval(isInject)<.05));
stack(2,2) = sum(pval(isInject)>.05);
stack(3,2) = sum((rho(isInject)>0) & (pval(isInject)<.05));
stack(:,2) = stack(:,2) / sum(stack(:,2));

barh(categorical(["Control","Inject"]), stack', 'stacked')
%ylim([0,sum(control_stack)])
xlabel('Proportion of clusters')

c = colororder;
textY = 2.5;
textX = stack(1,2)/2;
text(textX,textY, '-Corr','HorizontalAlignment','center','Color',c(1,:),'FontSize',20)
textX = stack(1,2) + stack(2,2)/2;
text(textX,textY, 'NS','HorizontalAlignment','center','Color',c(2,:),'FontSize',20)
textX = stack(1,2) + stack(2,2) + stack(3,2)/2;
text(textX,textY, '+Corr','HorizontalAlignment','center','Color',c(3,:),'FontSize',20)


%% Scatter plot correlation vs distance
figure(7); clf

subplot(2,1,1)
x = goodClusters.distFromInj(isControl);
y = rho(isControl);
scatter(x,y, 'filled', 'MarkerEdgeAlpha',0, 'MarkerFaceAlpha',.8)
xlim([600,1500])
title('Control')
ylabel('Correlation: Spiking vs concentration')

subplot(2,1,2)
x = goodClusters.distFromInj(isInject);
y = rho(isInject);
scatter(x,y, 'filled', 'MarkerEdgeAlpha',0, 'MarkerFaceAlpha',.8)
xlim([600,1500])
title('Inject')
ylabel('Correlation: spiking vs concentration')
xlabel('Distance from injection (um)')

%% Cross correlation
[r,lag] = xcorr(conc(1,:), spkZ(1,:));
r = nan(size(conc,1), length(r));

for i=1:size(conc,1)
    r(i,:) = xcorr(conc(i,:), spkZ(i,:), 'normalized');
end

%% Plot mean cross correlation (all)
figure(8); clf
x = lag*binWidth;
addShadedLine(x, r, 'k','All')
hold on
%ylim([0,1])
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
%ylim([0,1])
plot([0,0],ylim,'k--')
plot(xlim,[0,0], 'k--')
xlabel('Lag (s)')
ylabel('Correlation: spiking vs concentration')
legend("Control", "Inject")

%% Mountain plot of p-Value for specific lag
lagPoint = 90;
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