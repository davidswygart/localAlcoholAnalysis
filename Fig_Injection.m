%%
figFolder = 'C:\Users\david\Indiana University\Bryant, Kathleen - Acute Ethanol\Head-fixed mouse figures\1_Injection\linkedImages\';
figFolder = 'Downloads';
d_lowEk = load("diffusion_ek0p8.mat");
d_highEk = load("diffusion_ek2.mat");
load("goodClusters.mat")

controlColor = [65, 2, 87]/255;
injectColor = [217,95,2]/255;

PC1_color = [55,126,184]/255;
PC2_color = [228,26,28]/255;

positiveCorr_color = [216,179,101]/255;
negativeCorr_color = [90,180,172]/255;

injectBoxTransparency = 0.2;

isControl = contains(goodClusters.group,'control') | contains(goodClusters.group,'drink');
isInject = contains(goodClusters.group,'inject');
%% Bin spikes around injection 
binWidth=10;
%binEdges = (-60*10):binWidth:(60*12);
target = 'microInjectionStart';
binEdges = -600:binWidth:2220;
x_time = binEdges(1:end-1);
[spkCounts,  bpod]  = binAroundTarget(goodClusters, target, binEdges);
spkZ = zscore(spkCounts,0, 2);
spkZ = spkZ - mean(spkZ(:,x_time<0 & x_time>=-100), 2);

spkZ = spkZ(:, x_time>=-200 & x_time<=60*12);
x_time = x_time(x_time>=-200 & x_time<=60*12);
binEdges = [x_time, x_time(end)+binWidth];
%% spikes around injection (average)
figure(123); clf; hold on;

yLim = [-.75,.75];
ylim(yLim);
xlim([x_time(1),x_time(end)]);
r = rectangle(Position=[0,yLim(1),120,yLim(2)-yLim(1)], FaceColor=[0,0,0], EdgeColor='none', FaceAlpha=injectBoxTransparency)
% text(60,0.8,'Injection','HorizontalAlignment','center','VerticalAlignment','bottom')
addShadedLine(x_time,spkZ(isControl,:),{'Color', controlColor});
addShadedLine(x_time,spkZ(isInject,:),{'Color', injectColor});
yline(0)

xlabel('Time (s)')
ylabel('Firing (Z-score)')

% leg = legend('aCSF','aCSF + EtOH','Location','best');
% legend('boxoff')
% leg.ItemTokenSize = [5,4];

f = gcf;
f.Units = "inches";
f.Position = [2,2,          4,           1.8];
exportgraphics(gcf,[figFolder,'grandMean.pdf'], "ContentType","vector","BackgroundColor","none")

%% save for R
saveCsvForR_repeatedMeasuresMixedModel( ...
    {spkZ(isControl,:), spkZ(isInject,:)}, ...
    {"control","inject"}, ...
    'injectFiring.csv')


pVals = nan(size(spkZ,2),1);
for i=1:size(spkZ,2)
    [~,p]=ttest2(spkZ(isControl,i), spkZ(isInject,i),'Vartype','unequal');
    pVals(i) = p;
end
%h = fdr_bh(pVals, 0.05, 'dep');
%scatter(x_time(h), zeros(sum(h),1), "*")
%% Run PCA
weights = nan(size(spkZ,1),1);
weights(isControl) = 0.5 / sum(isControl);
weights(isInject) = 0.5 / sum(isInject);

% [coeff,score,latent,~,explained] = pca(spkZ'); % <-- Data are (time(bins) x neurons)
[coeff,score,latent,~,explained] = pca(spkZ', "VariableWeights",weights); % <-- Data are (time(bins) x neurons)
coeff = diag(sqrt(weights))*coeff;
%% Plot explained variance or Scree
figure(123);clf
% plot(latent,'.-')
ylabel('Eigenvalue')
bar(explained,'k');  
% ylabel(["Explained", "variance (%)"])
ylabel("Explained variance (%)")
xlabel('PC')
xlim([.50,3.5])
ylim([0,20])

figScale = 4;
a = gca;
a.FontSize = a.FontSize*figScale;
a.LineWidth = a.LineWidth*figScale;
f = gcf;
f.Units = "inches";
f.OuterPosition = [2,2,          .75*figScale,           1.5*figScale];
exportgraphics(gcf,[figFolder,'Scree.pdf'],"ContentType","vector",BackgroundColor="none")

%% Plot PC pattern (score)
figure(123);clf;hold on
% yLim = [-15,15];
yLim = [-.8,.8];
ylim(yLim);

plot(x_time, score(:,1),'LineWidth',1.5,'Color',PC1_color);
plot(x_time, score(:,2),'LineWidth',1.5,'Color',PC2_color);
% plot(x_time, score(:,3),'LineWidth',1.5,'Color',[0,.8,.2,.3]);

% plot([0,120], [14,14], 'Color',[.5,.5,.5], 'LineWidth',2)
rectangle(Position=[0,yLim(1),120,yLim(2)-yLim(1)], FaceColor=[0,0,0], EdgeColor='none',FaceAlpha=injectBoxTransparency)

xlabel('Time (s)'); ylabel('Firing (arbitrary units)');

hold on
yline(0)
%plotBpod(bpod)

% leg = legend({'PC1';'PC2'},'Location','best');
% legend('boxoff')
% leg.ItemTokenSize = [5,4];
xlim([x_time(1), x_time(end)])

f = gcf;
f.Units = "inches";
f.Position = [2,2,    4,         1.8];
exportgraphics(gcf,[figFolder,'pca_score.pdf'],"ContentType","vector",BackgroundColor="none")

%% PC1 loadings, split by group  (mountain plot)
figure(123);clf; hold on;
interestingPC = 1;

[f,x] = ecdf(coeff(isControl,interestingPC));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',1.5,'Color',controlColor)

[f,x] = ecdf(coeff(isInject,interestingPC));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',1.5,'Color',injectColor)

xline(0)
xlabel(['PC', num2str(interestingPC), ' loadings'])
ylabel('Probability')
% leg = legend("Control", "Inject","Location","northwest");
% legend('boxoff')
% leg.ItemTokenSize = [5,4];

xlim([-.1, 0.2])
ylim([0,.5])
box on

[h,p,ks2stat] = kstest2(coeff(isControl,interestingPC),coeff(isInject,interestingPC) );
% text(-.1,.4,['p=',num2str(p,2)], 'FontSize',5)
f = gcf;
f.Units = "inches";
f.Position = [2,2,          1.9,            1.5];
exportgraphics(gcf,[figFolder,'pca_PC1_loading.pdf'],"ContentType","vector",BackgroundColor="none")
%% PC2 loadings, split by group  (mountain plot)
figure(123);clf; hold on;
interestingPC = 2;

[f,x] = ecdf(coeff(isControl,interestingPC));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',1.5,'Color',controlColor)

[f,x] = ecdf(coeff(isInject,interestingPC));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',1.5,'Color',injectColor)

xline(0)
xlabel(['PC', num2str(interestingPC), ' loadings'])
ylabel('Probability')
% set(gca,"YTickLabel",['  '])

xlim([-.1, 0.17])
ylim([0,.5])
box on

[h,p,ks2stat] = kstest2(coeff(isControl,interestingPC),coeff(isInject,interestingPC) );
% text(-.1,.4,['p=',num2str(p,2)], 'FontSize',5)
f = gcf;
f.Units = "inches";
f.Position = [2,2,         1.9,            1.5];
exportgraphics(gcf,[figFolder,'pca_PC2_loading.pdf'],"ContentType","vector",BackgroundColor="none")

%% Combine high and low ek diffusion estimates

diffusion = struct;
diffusion.time = d_highEk.diffusion.time;
diffusion.dist = d_highEk.diffusion.dist;
diffusion.conc_mean = (d_highEk.diffusion.conc + d_lowEk.diffusion.conc)/2;
diffusion.conc_range = abs(d_highEk.diffusion.conc - d_lowEk.diffusion.conc)/2;


diffusion.time = [-200, diffusion.time];
diffusion.conc_mean = cat(2, zeros(size(diffusion.conc_mean, 1),1), diffusion.conc_mean);
diffusion.conc_range = cat(2, zeros(size(diffusion.conc_mean, 1),1), diffusion.conc_range);

%% Calculate the lowest 10%, median, and top 10% of cluster distances
sortedDist = sort(goodClusters.distFromInj(contains(goodClusters.group,'inject')));
sortedDist(round(length(sortedDist)*0.1));
sortedDist(round(length(sortedDist)*0.5));
sortedDist(round(length(sortedDist)*0.9));

%% plot concentration at different distances
figure(123); clf
hold on

yLim = [0,250];
ylim(yLim);
rectangle(Position=[0,yLim(1),120,yLim(2)-yLim(1)], FaceColor=[0,0,0], EdgeColor='none',FaceAlpha=injectBoxTransparency)

thisDist = 43;
shadedErrorBar(diffusion.time, diffusion.conc_mean(thisDist,:), diffusion.conc_range(thisDist,:), ...
    'LineProps', {'b-','LineWidth',1});

thisDist = 53;
shadedErrorBar(diffusion.time, diffusion.conc_mean(thisDist,:), diffusion.conc_range(thisDist,:), ...
    'LineProps', {'g-','LineWidth',1});

thisDist = 69;
shadedErrorBar(diffusion.time, diffusion.conc_mean(thisDist,:), diffusion.conc_range(thisDist,:), ...
    'LineProps', {'r-','LineWidth',1});


xlabel('Time (s)')
ylabel('Predicted [EtOH] (mg/dL)')
xlim([-200,720])

% leg = legend('430 um', '530 um', '690 um');
% legend('boxoff')
% leg.ItemTokenSize = [5,4];
f = gcf;
f.Units = "inches";
f.Position = [2,2,               2.5,            2];
exportgraphics(gcf,[figFolder,'concentration_atDistances.pdf'],"ContentType","vector","BackgroundColor","none")

%% plot histogram of cluster distances vs. peak concentration

figure(123);clf
hold on

% Plot peak ethanol
yyaxis left
[peakConc, inds] = max(diffusion.conc_mean,[],2);
linInd = sub2ind(size(diffusion.conc_mean), (1:length(inds))', inds);
peakConc_range =  diffusion.conc_range(linInd);
shadedErrorBar(diffusion.dist,peakConc,peakConc_range, ...
        'LineProps', {'-k','LineWidth',1});
ax = gca;
ax.YColor = 'k';
ylabel('Predicted peak [EtOH] (mg/dL)','Color', 'k')
% set(gca, 'YScale', 'log')
ylim([0,400])

% Plot Cluster distance
yyaxis right
b = 200:50:2000;
histogram(goodClusters.distFromInj(contains(goodClusters.group,'inject')),b)
% histogram(goodClusters.distFromInj,b)
ylabel('Number of neurons', 'Rotation',-90)
% set(get(gca,'YLabel'),'Rotation',-90)
% ylp = get(ylh, 'Position');
% ext=get(y_h,'Extent');
% set(y_h, 'Rotation',270, 'Position',ylp+[ext(3) 0 0])


% Set X axis
xlim([280,870])
xlabel('Distance from injection (µm)')


f = gcf;
f.Units = "inches";
f.Position = [2,2,2,2];
exportgraphics(gcf,[figFolder,'concentration_ClusterDistAndPeak.pdf'], "BackgroundColor","none","ContentType","vector")

%% Calculate predicted alcohol time per cluster
diffusion.conc = diffusion.conc_mean;
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
ylabel('Firing (Z-score)')

[p,S] = polyfit(x(:),y(:),2);
xp = 0:2:160;
[yp, delta] = polyval(p, xp, S);
hold on
shadedErrorBar(xp,yp,delta);
xlim([0,165])


%% Correlation of firing to concentration (lag = 0)
 [rho,pval] = corr(conc',spkZ');

 identity = logical(eye(length(rho)));
 pval = pval(identity);
 rho = rho(identity);

%% Corelation to firing (stacked bar) for paper
figure(123); clf; hold on; 
stack = nan(3,2);
stack(1,2) = sum(rho<0 & pval<.05 & isInject);
stack(2,2) = sum(pval>.05 & isInject);
stack(3,2) = sum(rho>0 & pval<.05 & isInject);


stack(1,1) = sum(rho<0 & pval<.05 & isControl);
stack(2,1) = sum(pval>.05 & isControl);
stack(3,1) = sum(rho>0 & pval<.05 & isControl);

p = twoProportionZtest(stack(:,1), stack(:,2));

stack(:,2) = stack(:,2) / sum(stack(:,2));
stack(:,1) = stack(:,1) / sum(stack(:,1));

% groupNames = ["aCSF","aCSF+EtOH"];
% x = categorical(groupNames);
% x = reordercats(x, groupNames);
b = bar([1,2], stack', 'stacked', '');
b(1).FaceColor = negativeCorr_color;
b(2).FaceColor = [.6,.6,.6];
b(3).FaceColor = positiveCorr_color;
%ylim([0,sum(control_stack)])
% set(gca, "XTickLabel", [])



% c = colororder;z
% textY = 2.7;
% textX = stack(1,2)/2;
% text(textX,textY, '-Corr','HorizontalAlignment','center','Color',c(1,:),'FontSize',8)
% textX = stack(1,2) + stack(2,2)/2;
% text(textX,textY, 'NS','HorizontalAlignment','center','Color',c(2,:),'FontSize',8)
% textX = stack(1,2) + stack(2,2) + stack(3,2)/2;
% text(textX,textY, '+Corr','HorizontalAlignment','center','Color',c(3,:),'FontSize',8)


ylabel('Proportion of neurons')

f = gcf;
%hA = axes(f);

box off
set(gca, "XTickLabel", [])
xlim([0.5,2.5])
% xtickangle(45)
% 
% leg = legend('-','NS', '+','Location','eastoutside');
% legend('boxoff')
% leg.ItemTokenSize = [5,4];


f.Units = "inches";
f.Position = [2,2,1.6,2.2];
exportgraphics(gcf,[figFolder,'concentration_stackedBar.pdf'],"ContentType","vector","BackgroundColor","none")

%% plot mean firing for high and low correlation clusters (Control)
figure(123); clf; hold on

yLim = [-1.1,2.7];
ylim(yLim);
rectangle(Position=[0,yLim(1),120,yLim(2)-yLim(1)], FaceColor=[0,0,0], EdgeColor='none',FaceAlpha=injectBoxTransparency)

y = spkZ(rho<0 & pval<.05 & isControl,:);
addShadedLine(x_time, y, {'Color', negativeCorr_color, 'Linewidth', 1});

y = spkZ(rho>0 & pval<.05 & isControl,:);
addShadedLine(x_time, y, {'Color', positiveCorr_color, 'Linewidth', 1});

xlim([x_time(1), x_time(end)])
yline(0)

% xlabel('Time (s)')
set(gca, "XTickLabels", [])
ylabel("Firing (Z-score)")
% title('Control')
% leg = legend('Control','Inject','Location','best');
% legend('boxoff')
% leg.ItemTokenSize = [5,4];

f = gcf;
f.Units = "inches";
f.Position = [2,2,            3,                 1.2];
exportgraphics(gcf,[figFolder,'Control_corr.pdf'], "ContentType","vector","BackgroundColor","none")

%% plot mean firing for high and low correlation clusters (Inject)
figure(123); clf; hold on

yLim = [-1.5,2.3];
ylim(yLim);
rectangle(Position=[0,yLim(1),120,yLim(2)-yLim(1)], FaceColor=[0,0,0], EdgeColor='none',FaceAlpha=injectBoxTransparency)

y = spkZ(rho<0 & pval<.05 & isInject,:);
addShadedLine(x_time, y,{'Color', negativeCorr_color, 'Linewidth', 1});

y = spkZ(rho>0 & pval<.05 & isInject,:);
addShadedLine(x_time, y,{'Color', positiveCorr_color, 'Linewidth', 1});

xlim([x_time(1), x_time(end)])
yline(0)
xlabel('Time (s)')
ylabel("Firing (Z-score)")
% title('Inject')
% leg = legend('Control','Inject','Location','best');
% legend('boxoff')
% leg.ItemTokenSize = [5,4];

f.Units = "inches";
f.Position = [2,2,             3,               1.4];
exportgraphics(gcf,[figFolder,'Inject_corr.pdf'], "ContentType","vector","BackgroundColor","none")
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
[r_temp,lag] = xcorr(spkZ(1,:), conc(1,:)); % This is just to get lag vector and size of correlation vector
lag = lag*binWidth;
r = nan(size(conc,1), length(r_temp));  % empty matrix to store correlation values
r_shuff = r;  % empty matrix to store shuffled correlation values
r_auto_conc = r;  % empty matrix to store [EtOH] autocorrelation
r_auto_spike = r;  % empty matrix to store spike autocorrelation

for i=1:size(conc,1) % loop through every neuron to get xcorr between spikes and predicted concentration
    x = spkZ(i,:) - mean(spkZ(i,:));
    y = conc(i,:);
    r(i,:) = xcorr(x,y, 'normalized'); % get real cross-correlation
    r_auto_conc(i,:) = xcorr(y,y, 'normalized'); % 
    r_auto_spike(i,:) = xcorr(x,x, 'normalized'); % 

    x = x(randperm(length(x))); % shuffle spikes for control
    r_shuff(i,:) = xcorr(x,y,'normalized'); % get shuffled cross-correlation
end

%r_auto_conc = r_auto_conc - r; % corrected autocorrelation
%r_auto_spike = r_auto_spike - r; % corrected autocorrelation

%% Plot mean cross correlation (split by type)
figure(123); clf
hold on

% plot shuffled controls 
addShadedLine(lag, r_shuff(isControl,:), {'--', 'Color', (1-controlColor)*.3 + controlColor});
addShadedLine(lag, r_shuff(isInject,:), {'--', 'Color', (1-injectColor)*.3 + injectColor});

% plot real data
addShadedLine(lag, r(isControl,:), {'Color', controlColor});
addShadedLine(lag, r(isInject,:), {'Color', injectColor});

yline(0)
xline(0)

xlabel('Lag (s)')
ylabel('Spiking correlation to [EtOH]')
xlim([-500,700])

f = gcf;
f.Units = "inches";
f.Position = [2,2,            3,             2.4];
exportgraphics(gcf,[figFolder,'crossCorrelation.pdf'], "ContentType","vector","BackgroundColor","none")

%% Plot mean cross correlation (EtOH only) compared to corrected autocorrelation
fig =figure(123); clf
set(fig,'defaultLegendAutoUpdate','off');
hold on

% plot shuffled controls 
addShadedLine(lag, r_shuff(isInject,:), {'--', 'Color', (1-injectColor)*.3 + injectColor});

% plot corrected autocorrelation  
addShadedLine(lag, r_auto_conc(isInject,:), {'-', 'Color', [.9,.9,.5]});
addShadedLine(lag, r_auto_spike(isInject,:), {'-', 'Color',[.9,.1,.9]});

% plot real data
addShadedLine(lag, r(isInject,:), {'Color', injectColor});

legend('Shuffled', 'Auto ([EtOH])', 'Auto (spikes)', 'Cross correlation')


yline(0)
xline(0)

xlabel('Lag (s)')
ylabel('Spiking correlation to [EtOH]')
xlim([-500,700])

f = gcf;
f.Units = "inches";
%f.Position = [2,2,            3,             2.4];
exportgraphics(gcf,[figFolder,'crossCorrelation.pdf'], "ContentType","vector","BackgroundColor","none")
%% Mountain plot of correlation for peak min lag
[~, lagInd] = min(mean(r(isInject,:)));

figure(123); clf; hold on;
[f,x] = ecdf(r(isControl,lagInd));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',1.5,'Color',controlColor)
hold on

[f,x] = ecdf(r(isInject,lagInd));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',1.5,'Color',injectColor)

xlabel('Correlation')
ylabel('Probability')

% title(['Spike lag ', num2str(lag(lagInd)), 's'])
xline(0)

xlim([-.6, 0.5])
ylim([0,.5])
box on

f = gcf;
f.Units = "inches";
f.Position = [2,2,          1.6,            1.1];
exportgraphics(gcf,[figFolder,'correlation_peakMinLag.pdf'], "ContentType","vector","BackgroundColor","none")


[~,pval,~,ttestStats] = ttest2(r(isControl,lagInd), r(isInject,lagInd), "Vartype","unequal")

%% Mountain plot of correlation for peak max lag
[~, lagInd] = max(mean(r(isInject,:)));

figure(123); clf; hold on;
xline(0)

[f,x] = ecdf(r(isControl,lagInd));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',1.5,'Color',controlColor)

[f,x] = ecdf(r(isInject,lagInd));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',1.5,'Color',injectColor)

xlabel('Correlation')
ylabel('Probability')

xlim([-.4, 0.7])
ylim([0,.5])
box on

f = gcf;
f.Units = "inches";
f.Position = [2,2,          1.6,            1.1];
exportgraphics(gcf,[figFolder,'correlation_peakMaxLag.pdf'], "ContentType","vector","BackgroundColor","none")

[~,pval,~,ttestStats] = ttest2(r(isControl,lagInd), r(isInject,lagInd), "Vartype","unequal")
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