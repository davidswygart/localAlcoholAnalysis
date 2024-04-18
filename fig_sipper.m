figFolder = 'C:\Users\david\OneDrive - Indiana University\localAlcohol\Figures\2_Sipper\matlabExports\';
% figFolder = 'C:\Users\dis006\OneDrive - Indiana University\localAlcohol\Figures\2_Sipper\matlabExports\';
load("goodClusters.mat")
load("group.mat")

controlColor = [65, 2, 87]/255;
drinkColor = [5, 105, 67]/255;
injectColor = [217,95,2]/255;

PC1_color = [55,126,184]/255;
PC2_color = [228,26,28]/255;

isControl = contains(goodClusters.group,'control');
isDrink = contains(goodClusters.group,'drink');
isInject = contains(goodClusters.group,'inject');


injectBoxTransparency = 0.2;
%% Fluid consumed
group.mlPerkg = 1000 * group.fluidConsumed ./ group.mouseWeight;

controlConsumed = group.mlPerkg(contains(group.group,'control'));
injectedConsumed = group.mlPerkg(contains(group.group, 'injected'));
drinkConsumed = group.mlPerkg(contains(group.group, 'drink'));


figure(123);clf; hold on
groupNames = ["Control", "EtOH inject","EtOH sipper"];
x = categorical(groupNames);
x = reordercats(x,groupNames);
y = [mean(controlConsumed), mean(injectedConsumed),mean(drinkConsumed)];
err = [std(controlConsumed)/sqrt(length(controlConsumed)), std(injectedConsumed)/sqrt(length(injectedConsumed)), std(drinkConsumed)/sqrt(length(drinkConsumed))];

b = bar(x,y);
b.FaceColor = 'flat';
colorMult = 1.2;
b.CData(1,:) = controlColor;
b.CData(2,:) = injectColor;
b.CData(3,:) = drinkColor;
errorbar(x,y,err,'k', 'LineStyle','none')
% 
scatter(ones(length(controlConsumed),1)*.75 + .5*rand(length(controlConsumed),1) , controlConsumed, 'k.')
scatter(ones(length(injectedConsumed),1)*1.75 + .5*rand(length(injectedConsumed),1) , injectedConsumed, 'k.')
scatter(ones(length(drinkConsumed),1)*2.75 + .5*rand(length(drinkConsumed),1) , drinkConsumed, 'k.')

ylabel('Fluid consumed (ml/kg)')

f = gcf;
f.Units = "inches";
f.Position = [2,2,           1.2,                1.7];
exportgraphics(gcf,[figFolder,'volumeConsumed.pdf'],"BackgroundColor","none","ContentType","vector")


%% Bin spikes around drinking  Full time;
binWidth=10;
target = 'sipperStart';
binEdges = -1320:binWidth:1500;
x_fullTime = binEdges(1:end-1);

[spkCounts,  bpod]  = binAroundTarget(goodClusters, target, binEdges);%, 'smooth');
spkZ_fullTime = zscore(spkCounts,0, 2);
spkZ_fullTime = spkZ_fullTime - mean(spkZ_fullTime(:,x_fullTime<0 & x_fullTime>=-100), 2);

startTime = -200;
stopTime = 1400;
spkZ_fullTime = spkZ_fullTime(:, x_fullTime>=startTime & x_fullTime<=stopTime);
x_fullTime = x_fullTime(x_fullTime>=startTime & x_fullTime<=stopTime);

%% spikes around drinking Full time (average)
figure(123); clf; hold on
yLim = [-.5,1];
ylim(yLim);

xlim([x_fullTime(1), x_fullTime(end)])

addShadedLine(x_fullTime,spkZ_fullTime(isControl,:),{'Color', controlColor});
addShadedLine(x_fullTime,spkZ_fullTime(isInject,:),{'Color', injectColor});
addShadedLine(x_fullTime,spkZ_fullTime(isDrink,:),{'Color', drinkColor});
yline(0)

% plot([0,900], [1,1], 'Color',[.5,.5,.5], 'LineWidth',2)
rectangle(Position=[0,yLim(1),900,yLim(2)-yLim(1)], FaceColor=[0,0,0,injectBoxTransparency], EdgeColor='none')
text(450,1,"Sipper active","HorizontalAlignment","center","VerticalAlignment","top")

xlabel('Time (s)')
ylabel('Firing (Z-score)')

% leg = legend('Control','EtOH consumed','EtOH injected','Location','northeast');
% legend('boxoff')
% leg.ItemTokenSize = [5,4];

f = gcf;
f.Units = "inches";
f.Position = [2,2,2.8,1.8];
exportgraphics(gcf,[figFolder,'grandMean.pdf'], "ContentType","vector","BackgroundColor","none")

%% Spikes around sipper valve
binWidth=0.1;
binEdges = -1320:binWidth:1500;
target = 'sipperStart';
[tempCounts,  ~]  = binAroundTarget(goodClusters, target, binEdges);
clusterSTD = std(tempCounts, 0 ,2);
clusterMean = mean(tempCounts,2);

binEdges = -2:binWidth:14;
x_time = binEdges(1:end-1);
spkCounts  = binAroundValve(goodClusters, binEdges);%,'smooth');
spkAvg = mean(spkCounts ,3);
spkZ = (spkAvg-clusterMean) ./  clusterSTD;
spkZ = spkZ - mean(spkZ(:, x_time<0 & x_time>=-1),2);
%% Spikes around sipper valve - mean
figure(123); clf; hold on

ylim([-.15,.35])
plot([0,0],ylim,'--','Color',[.5,.5,.5])
plot([10.1,10.1],ylim,'--','Color',[.5,.5,.5])

addShadedLine(x_time, spkZ(isControl,:), {'Color', controlColor});
addShadedLine(x_time, spkZ(isInject,:), {'Color', injectColor});
addShadedLine(x_time, spkZ(isDrink,:), {'Color', drinkColor});
yline(0)

ylabel('Firing (Z-score)')
% xlabel('Time (s)')
xlabel(' ')
set(gca, "XTickLabel", [])
xlim([x_time(1), x_time(end)])

% leg = legend('Control','EtOH consumed','EtOH injected',Location='north');
% legend('boxoff')
% leg.ItemTokenSize = [5,4];

f = gcf;
f.Units = "inches";
f.Position = [2,2,           4.2,               1.8];
exportgraphics(gcf,[figFolder,'valveMean.pdf'],"ContentType","vector","BackgroundColor","none")
%% Run PCA
weights = nan(size(spkZ,1),1);
weights(isControl) = 1/3 / sum(isControl);
weights(isInject) = 1/3 / sum(isInject);
weights(isDrink) = 1/3 / sum(isDrink);

% [coeff,score,latent,~,explained] = pca(spkZ'); % <-- Data are (time(bins) x neurons)
[coeff,score,latent,~,explained] = pca(spkZ', "VariableWeights",weights); % <-- Data are (time(bins) x neurons)
coeff = diag(sqrt(weights))*coeff;
%% Plot explained variance or Scree
figure(123);clf
% plot(latent,'.-')
bar(explained,'k');  
ylabel(["Explained", "variance (%)"])
xlabel('PC')
xlim([.50,3.5])

f = gcf;
f.Units = "inches";
f.Position = [2,2, 0.9, .9];
exportgraphics(gcf,[figFolder,'pca_scree.pdf'],"ContentType","vector","BackgroundColor","none")

%% Plot PC pattern (score)
figure(123); clf
hold on

% ylim([-3,8])
ylim([-.17,.38])
plot([0,0],ylim,'--','Color',[.5,.5,.5])
plot([10.1,10.1],ylim,'--','Color',[.5,.5,.5])

plot(x_time, score(:,1),'LineWidth',1,'Color',PC1_color);
plot(x_time, score(:,2),'LineWidth',1,'Color',PC2_color);
% plot(x_time, score(:,3),'LineWidth',2,'Color',[0,.8,.2,.5]);
%plot(x_time, score(:,4),'LineWidth',2,'Color',[.5,.5,.5,.5]);

xlabel('Time (s)');
ylabel('Firing (arbitrary units)');
xlim([x_time(1), x_time(end)])

hold on
yline(0)

% leg = legend({'PC1';'PC2'},'Location','north');
% legend('boxoff')
% leg.ItemTokenSize = [5,4];


f = gcf;
f.Units = "inches";
f.Position = [2,2,        4.2,             1.8];
exportgraphics(gcf,[figFolder,'pca_score.pdf'], "ContentType","vector", "BackgroundColor","none")

%% PC1 loadings, split by group  (mountain plot)
figure(123); clf; hold on;
interestingPC = 1;

xline(0)

[f,x] = ecdf(coeff(isControl,interestingPC));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',1.5,'Color',controlColor)

[f,x] = ecdf(coeff(isInject,interestingPC));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',1.5,'Color', injectColor)

[f,x] = ecdf(coeff(isDrink,interestingPC));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',1.5,'Color',drinkColor)

xlabel(['PC', num2str(interestingPC), ' loadings'])
xlim([-.05,0.25])
ylabel('Probability')

[h,p] = kstest2(coeff(isControl,interestingPC),coeff(isDrink,interestingPC) );
% text(-.14,.4,['CvD p=',num2str(p*2,2)], 'FontSize',5)

[h,p] = kstest2(coeff(isControl,interestingPC),coeff(isInject,interestingPC) );
% text(-.14,.35,['CvI p=',num2str(p*2,2)], 'FontSize',5)

% pc1_Thresh = -.05;
% plot([pc1_Thresh,pc1_Thresh], [0,.5], '-.k')
% leg = legend('Control','EtOH consumed','EtOH injected',"Location","northeast");
% legend('boxoff')
% leg.ItemTokenSize = [5,4];\

box on

f = gcf;
f.Units = "inches";
f.Position = [2,2,2,1.2];
exportgraphics(gcf,[figFolder,'pca_PC1Mountain.pdf'])
%% PC2 loadings, split by group  (mountain plot)
figure(123); clf; hold on;
interestingPC = 2;

xline(0)

[f,x] = ecdf(coeff(isControl,interestingPC));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',1.5,'Color',controlColor)

[f,x] = ecdf(coeff(isInject,interestingPC));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',1.5,'Color', injectColor)

[f,x] = ecdf(coeff(isDrink,interestingPC));
f(f>0.5) = 1 - f(f>0.5);
plot(x,f,'LineWidth',1.5,'Color',drinkColor)

xlabel(['PC', num2str(interestingPC), ' loadings'])
xlim([-.1,.2])
% ylabel('Folded probability')

[h,p] = kstest2(coeff(isControl,interestingPC),coeff(isDrink,interestingPC) );
% text(-.14,.4,['CvD p=',num2str(p*2,2)], 'FontSize',5)
set(gca,'Yticklabel',[]) 

[h,p] = kstest2(coeff(isControl,interestingPC),coeff(isInject,interestingPC) );
% text(-.14,.35,['CvI p=',num2str(p*2,2)], 'FontSize',5)

% pc1_Thresh = -.05;
% plot([pc1_Thresh,pc1_Thresh], [0,.5], '-.k')
% leg = legend('Control','EtOH consumed','EtOH injected',"Location","northwest");
% legend('boxoff')
% leg.ItemTokenSize = [5,4];

box on

f = gcf;
f.Units = "inches";
f.Position = [2,2,2,1.2];
exportgraphics(gcf,[figFolder,'pca_PC2Mountain.pdf'])

%% Plot PC1 vs. PC2
figure(11)
scatter(coeff(isInject,1), coeff(isInject,2), 'filled')
yline(0)
xline(0)
title('Inject')

figure(12)
scatter(coeff(isDrink,1), coeff(isDrink,2), 'filled')
yline(0)
xline(0)
title('Drink')

figure(13)
scatter(coeff(isControl,1), coeff(isControl,2), 'filled')
yline(0)
xline(0)
title('Control')
