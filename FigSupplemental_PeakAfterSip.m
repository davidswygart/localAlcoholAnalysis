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

%% Spikes around sipper valve - Group by N valve
binWidth=0.1;
binEdges = -1320:binWidth:1500;
target = 'sipperStart';
[tempCounts,  ~]  = binAroundTarget(goodClusters, target, binEdges);
clusterSTD = std(tempCounts, 0 ,2);
clusterMean = mean(tempCounts,2);

binEdges = -2:binWidth:14;
x_time = binEdges(1:end-1);
spkCounts  = binAroundValve(goodClusters, binEdges);%,'smooth');

valveBinning = 12;
clumped = nan(size(spkCounts,1),size(spkCounts,2),size(spkCounts,3)/valveBinning);
for i=1:size(spkCounts,3)/valveBinning
    ind1 = (i-1)*valveBinning+1;
    ind2 = ind1 + valveBinning - 1;
    clumped(:,:,i) = mean(spkCounts(:,:,ind1:ind2), 3);
end

spkZ_clumped = zscore(clumped,0,2);

spkZ = (clumped-clusterMean) ./  clusterSTD;
spkZ = spkZ - mean(spkZ(:, x_time<0 & x_time>=-1),2);

%% Plot traces for each valve opening with clumped valves
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
inhWindow = val2Ind(x_time, [1.5,3.4]); %time below 0
inh = squeeze(mean(spkZ_clumped(:,inhWindow(1):inhWindow(2),:), 2));


figure(7);clf

hold on
title('"Inh" peak (1.6-3.2s after valve)')
addShadedLine([],inh(isControl,:), {'Color', controlColor});
addShadedLine([],inh(isDrink,:),{'Color', drinkColor});
addShadedLine([],inh(isInject,:),{'Color', injectColor});
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
addShadedLine(x_time, y, {'Color',controlColor});
y = spkZ_clumped(isDrink,:,1);
addShadedLine(x_time, y, {'Color', drinkColor});
y = spkZ_clumped(isInject,:,1);
addShadedLine(x_time, y, {'Color', injectColor});

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
addShadedLine(x_time, y, {'Color',controlColor});
y = spkZ_clumped(isDrink,:,3);
addShadedLine(x_time, y, {'Color',drinkColor});
y = spkZ_clumped(isInject,:,3);
addShadedLine(x_time, y, {'Color',injectColor});

yline(0)
plot([0,0],ylim,'--','Color',[.5,.5,.5])
plot([10.1,10.1],ylim,'--','Color',[.5,.5,.5])
xlim([x_time(1), x_time(end)])
ylabel('zscore')

nexttile
hold on
y = spkZ_clumped(isControl,:,5);
addShadedLine(x_time, y,  {'Color',controlColor});
y = spkZ_clumped(isDrink,:,5);
addShadedLine(x_time, y, {'Color',drinkColor});
y = spkZ_clumped(isInject,:,5);
addShadedLine(x_time, y, {'Color',injectColor});

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


%%
function vals = val2Ind(x, vals)
    for i=1:length(vals)
        [~,minInd]=min(abs(x-vals(i)));
        vals(i)=minInd;
    end
end
