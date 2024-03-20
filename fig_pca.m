%% spikes around injection 
binWidth=10;
binEdges = (-60*10):binWidth:(60*12);
target = 'microInjectionStart';
cLabel = 'zscore';
cRange = [-.4,6];

[spkCounts,  bpod]  = binAroundTarget(goodClusters, target, binEdges, 'smooth');
spkZ = zscore(spkCounts,0, 2);

figure(1); clf
plotSpikeHeatmapWithBpod(spkZ, binEdges, bpod, cLabel,cRange)
xlabel('Time (10s bins)')

%% Run PCA
[coeff,score,latent,~,explained] = pca(spkZ'); % <-- Data are (time(bins) x neurons)

%% Plot explained variance or Scree
figure(2);clf
plot(latent,'.-')
ylabel('Eigenvalue')
bar(explained,'k');  
% bar(cumsum(explained),'k');
%ylabel('Explained variance (%)')
xlabel('PC')

%% Plot PC pattern (score)
figure(3); clf

pcs = 1:3;
plot(binEdges(1:end-1), score(:,pcs));
xlabel('Time (10s bins)'); ylabel('Z-score');

hold on
plotBpod(bpod)

legend({'PC1';'PC2';'PC3'})
xlim([binEdges(1), binEdges(end)])
%% Identify groups
isControl = contains(goodClusters.group,'control') | contains(goodClusters.group,'drink');
isInject = contains(goodClusters.group,'inject');
%% PC2 loadings, split by group (histogram)
loadingThresh = 0.01;

figure(4)
clf
coeff_binEdges = -.15:.01:.15;
x = coeff_binEdges(1:end-1) + (coeff_binEdges(2)-coeff_binEdges(1))/2;
group = isControl;
control = histcounts(coeff(group,2), coeff_binEdges);
b = bar(x, control/sum(control), 'blue');
b.FaceAlpha = .5;
hold on
group = isInject;
injected = histcounts(coeff(group,2), coeff_binEdges);
b = bar(x, injected/sum(injected), 'red');
b.FaceAlpha = .5;

% % plot threshold
% plot([-loadingThresh,-loadingThresh], ylim, 'k--')
% plot([loadingThresh,loadingThresh], ylim, 'k--')

legend('control','injected')
xlabel('PC2 loading (coeff)')
ylabel('Count (norm)')

%% PC2 loadings, split by group  (mountain plot)
figure(4); clf
[f,x] = ecdf(coeff(isControl,2));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f)
hold on
[f,x] = ecdf(coeff(isInject,2));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f)

plot([0,0],[0,0.5], '--k')
xlabel('PC2 loading (coeff)')
ylabel('Folded probability')
legend("Control", "Inject","Location","northwest")

[h,p] = kstest2(coeff(isControl,2),coeff(isInject,2) );
text(-.12,.25,['p=',num2str(p,2)])

%% Plot spikes sorted by PC2
[sortedCoeff,sortInd] = sort(coeff(:,2));
spkZ_sorted = spkZ(sortInd,:);

figure(5); clf
plotSpikeHeatmapWithBpod(spkZ_sorted, binEdges, bpod, cLabel,cRange)
title('Sorted by PC2')
xlabel('Time (10s bins)')


%plot threshold limits
hold on
[~,threshInd]=min(abs(-sortedCoeff-loadingThresh));
plot(xlim,[threshInd,threshInd], 'r--')
[~,threshInd]=min(abs(sortedCoeff-loadingThresh));
plot(xlim,[threshInd,threshInd], 'r--')

%% Plot stacked bar graphs of loader N
controlLoaders = nan(1,3);
controlLoaders(1) = sum(coeff(isControl,2) < -loadingThresh);
controlLoaders(3) = sum(coeff(isControl,2) > loadingThresh);
controlLoaders(2) = sum(isControl) - controlLoaders(3) - controlLoaders(1);

injectLoaders = nan(1,3);
injectLoaders(1) = sum(coeff(isInject,2) < -loadingThresh);
injectLoaders(3) = sum(coeff(isInject,2) > loadingThresh);
injectLoaders(2) = sum(isInject) - injectLoaders(3) - injectLoaders(1);

figure(6); clf
subplot(1,2,1)
bar("Control",controlLoaders, 'stacked')
ylabel('Count')

subplot(1,2,2)
bar("Inject", injectLoaders, 'stacked')
legend('negative', 'none' , 'positive', 'Location','northeast')
%% plot mean split out by loaders
pLoad_control = spkZ(isControl & coeff(:,2)>loadingThresh, :);
pLoad_inject = spkZ(isInject & coeff(:,2)>loadingThresh, :);

nLoad_control = spkZ(isControl & coeff(:,2)<-loadingThresh, :);
nLoad_inject = spkZ(isInject & coeff(:,2)<-loadingThresh, :);

figure(7); clf

subplot(2,1,1)
x_time = binEdges(1:end-1) + (binEdges(2)-binEdges(1))/2;
addShadedLine(x_time,pLoad_control,'b','Control')
hold on
addShadedLine(x_time,pLoad_inject,'r','Inject')
plotBpod(bpod)
title('Positive loaders')
xlim([x_time(1), x_time(end)])

subplot(2,1,2)
addShadedLine(x_time,nLoad_control,'b','Control')
hold on
addShadedLine(x_time,nLoad_inject,'r','Inject')
plotBpod(bpod)
title('Negative loaders')
xlim([x_time(1), x_time(end)])
legend("Control","Inject")


%% %%%%%%%%%%%%%% alternative plots %%%%%%%%%%%%%%
%% Plot sorted coeff
figure(6); clf
plot(sortedCoeff)
hold on
plot(xlim,[0,0],'k--')
xlabel('sorted cluster #')
ylabel('PC2 loading (coeff)')

%plot threshold limits
hold on
[~,threshInd]=min(abs(-sortedCoeff-loadingThresh));
plot([threshInd,threshInd],ylim, 'r--')
[~,threshInd]=min(abs(sortedCoeff-loadingThresh));
plot([threshInd,threshInd],ylim, 'r--')


%% PCA loadings, split by group (3D)
figure(4);clf
group = isControl;
scatter3(coeff(group,1),coeff(group,2),coeff(group,3), 'b.')
xlabel('pc1')
ylabel('pc2')
zlabel('pc3')
hold on
group = isInject;
scatter3(coeff(group,1),coeff(group,2),coeff(group,3), 'r.')

%% PCA loadings, split by group (2D)
figure(4)
clf
subplot(1,3,1)
group = isControl;
scatter(coeff(group,1),coeff(group,2), 'b.')
hold on
group = isInject;
scatter(coeff(group,1),coeff(group,2), 'r.')
xlabel('pc1')
ylabel('pc2')

subplot(1,3,2)
group = isControl;
scatter(coeff(group,1),coeff(group,3), 'b.')
hold on
group = isInject;
scatter(coeff(group,1),coeff(group,3), 'r.')
xlabel('pc1')
ylabel('pc3')

subplot(1,3,3)
group = isControl;
scatter(coeff(group,2),coeff(group,3), 'b.')
hold on
group = isInject;
scatter(coeff(group,2),coeff(group,3), 'r.')
xlabel('pc2')
ylabel('pc3')

legend('Control','Inject')

%%
%%
%% Examine vars created by PCA
pc = 1; 

%% coeff is the weight for each neuron for each PC
figure(2)
imagesc(coeff);     
%%
plot(coeff(:,pc));              % The coefficients for PC1 
%%
plot(sort(coeff(:,pc)),'ko');   % Easier to make sense of sorted coefficients

%%
% split pos and neg loaders, these have different, often mirror-like, patterns
k = find(coeff(:,pc)>0);
plot(spkZ(:,k)*coeff(k,pc),'b');hold on;

k = find(coeff(:,pc)<0);
plot(spkZ(:,k)*coeff(k,pc),'r');hold on;
legend('positive loader','negative loader')

%% Can sort the matrix by the coeff for a PC to visualize how the PCA 
pc = 2
figure(3)
clim=[-3 3];
[k,l] = sort(coeff(:,pc));
subplot(2,3,1)
bar(explained(1:10));
xlabel('PC Number'); ylabel('Exp variance (%)');
subplot(2,3,4)
plot(score(:,1:3));
xlabel('Time (bins)'); ylabel('Ca activity');
legend({'PC1';'PC2';'PC3'})
axis tight
subplot(2,3,[2 5])
imagesc(zscore(spkZ)',clim);
xlabel('Time (bins)')
ylabel('Neuron')
title('Unsorted data');
subplot(2,3,[3 6])
imagesc(zscore(spkZ(:,l))',clim);
xlabel('Time (bins)')
ylabel('Neuron')
title('Sorted data');


