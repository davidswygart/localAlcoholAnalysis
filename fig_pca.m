%% spikes around injection 
binWidth=10;
binEdges = (-60*10):binWidth:(60*12);
target = 'microInjectionStart';
cLabel = 'zscore';
cRange = [-.5,8];

[spkCounts,  bpod]  = binAroundTarget(goodClusters, target, binEdges, 'smooth');
spkZ = zscore(spkCounts,0, 2);

clf; figure(1)
plotSpikeHeatmapWithBpod(spkZ, binEdges, bpod, cLabel,cRange)
title('Inject')
xlabel('Time (10s bins)')

%% Run PCA

[coeff,score,latent,~,explained] = pca(spkZ'); % <-- Data are (time(bins) x neurons)
spkZ = spkZ';
%% Examine vars created by PCA
pc = 1; % Looking at PC1

% coeff is the weight for each neuron for each PC
imagesc(coeff);               
plot(coeff(:,pc));              % The coefficients for PC1 
plot(sort(coeff(:,pc)),'ko');   % Easier to make sense of sorted coefficients

% project data on PC1
plot(spkZ*coeff(:,pc),'k');hold on; % talk about .* (element-wise) vs. * (matrix)

% split pos and neg loaders, these have different, often mirror-like, patterns
k = find(coeff(:,pc)>0);
plot(spkZ(:,k)*coeff(k,pc),'b');hold on;
coeff(k,pc)
k = find(coeff(:,pc)<0);
plot(spkZ(:,k)*coeff(k,pc),'r');hold on;
coeff(k,pc)

% score is the PC
plot(score(:,[1:3]));hold on
plot(spkZ*coeff(:,pc),'k--');hold on;

% latent is the eigenvalue (some folks argue eigs < 1 are useless, not me)
plot(latent);                       
plot((latent./sum(latent)).*100, 'ko');hold on;

% explained is the explained variance of each PC
bar(explained,'k');               
bar(cumsum(explained),'k');    % common to show exp variance like this  

%% Can sort the matrix by the coeff for a PC to visualize how the PCA 

figure()
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

%%
isControl = contains(goodClusters.group,'control') | contains(goodClusters.group,'drink');
control = spkZ(isControl,:);
isInject = contains(goodClusters.group,'inject');
inject = spkZ(isInject,:);