%load("allClusters.mat")

allClusters = allClusters(allClusters.fr>0.1, :);
allClusters = allClusters(~contains(allClusters.phyLabel,'noise'), :);
allClusters = allClusters(allClusters.presenceRatio>0.9, :);
allClusters = allClusters(allClusters.isiViolations<1, :);

%% choose bin width
binWidth=0.5;
binEdges = -3:binWidth:15;

%%
includedValves = 1:60;
spkCounts  = binAroundValve(allClusters, binEdges);%,'smooth');
spkAvg = mean(spkCounts(:,:,includedValves) ,3);

spkRate = spkAvg / binWidth;
spkZ = zscore(spkAvg,0, 2);
%spkZ = spkRate;

yname = 'zscore';
control = spkZ(contains(allClusters.group,'control'), :);
drink  = spkZ(contains(allClusters.group,'drink'), :);
inject = spkZ(contains(allClusters.group,'injected'), :);


%%
avgControl = mean(control,1);
stdControl = std(control,0, 1);
semControl = stdControl / sqrt(size(control,1));

avgDrink = mean(drink,1);
stdDrink = std(drink,0, 1);
semDrink = stdDrink / sqrt(size(drink,1));

avgInject = mean(inject,1);
stdInject = std(inject,0, 1);
semInject = stdInject / sqrt(size(inject,1));

figure(4)
clf
title('Average spiking after valve open')
shadedErrorBar(binEdges(1:end-1), avgControl,semControl, 'lineProps', 'b')
hold on
shadedErrorBar(binEdges(1:end-1), avgDrink,semDrink, 'lineProps', 'r')
shadedErrorBar(binEdges(1:end-1), avgInject,semInject, 'lineProps', 'g')

legend('control','drink','inject')
ylabel(yname)
xlabel('time (s)')
xlim([binEdges(1), binEdges(end)])
hold off



%%

[~,indPeakExc] = min(abs(binEdges-0.5));
[~,indPeakInh] = min(abs(binEdges-2));

valveBinning = 10;
clumpedSpk = nan(size(spkCounts,1),size(spkCounts,2),size(spkCounts,3)/valveBinning);
for i=1:size(spkCounts,3)/valveBinning
    ind1 = (i-1)*valveBinning+1;
    ind2 = ind1 + valveBinning - 1;
    clumpedSpk(:,:,i) = mean(spkCounts(:,:,ind1:ind2), 3);

end

%%
%spkZ = zscore(spkCounts,0, 2);
spkZ = zscore(clumpedSpk,0, 2);



controlExc = spkZ(contains(allClusters.group,'control'), indPeakExc, :);
drinkExc  = spkZ(contains(allClusters.group,'drink'), indPeakExc, :);
injectExc = spkZ(contains(allClusters.group,'injected'), indPeakExc, :);

controlInh = spkZ(contains(allClusters.group,'control'), indPeakInh, :);
drinkInh  = spkZ(contains(allClusters.group,'drink'), indPeakInh, :);
injectInh = spkZ(contains(allClusters.group,'injected'), indPeakInh, :);

avgControlExc = squeeze(mean(controlExc,1));
avgDrinkExc = squeeze(mean(drinkExc,1));
avgInjectExc = squeeze(mean(injectExc,1));

semControlExc = squeeze(std(controlExc,0,1)) ./ sqrt(size(controlExc,1));
semDrinkExc = squeeze(std(drinkExc,0,1)) ./ sqrt(size(drinkExc,1));
semInjectExc = squeeze(std(injectExc,0,1)) ./ sqrt(size(injectExc,1));

avgControlInh = squeeze(mean(controlInh,1));
avgDrinkInh = squeeze(mean(drinkInh,1));
avgInjectInh = squeeze(mean(injectInh,1));

semControlInh = squeeze(std(controlInh,0,1)) ./ sqrt(size(controlInh,1));
semDrinkInh = squeeze(std(drinkInh,0,1)) ./ sqrt(size(drinkInh,1));
semInjectInh = squeeze(std(injectInh,0,1)) ./ sqrt(size(injectInh,1));

%time = linspace(0,length(includedValves)*15.1,length(includedValves));
time = includedValves(1:valveBinning:end);

figure(5)
clf

subplot(2,1,1)
hold on
title('"Exc" peak (0.5-1s after valve)')
shadedErrorBar(time, avgControlExc, semControlExc, 'lineProps', 'b')
shadedErrorBar(time, avgDrinkExc, semDrinkExc, 'lineProps', 'r')
shadedErrorBar(time, avgInjectExc, semInjectExc, 'lineProps', 'g')
plot([time(1),time(end)], [0,0], '--k')
ylabel('zscore')
hold off

subplot(2,1,2)
hold on
title('"Inh" peak (2-2.5s after valve)')
shadedErrorBar(time, avgControlInh, semControlInh, 'lineProps', 'b')
shadedErrorBar(time, avgDrinkInh, semDrinkInh, 'lineProps', 'r')
shadedErrorBar(time, avgInjectInh, semInjectInh, 'lineProps', 'g')
plot([time(1),time(end)], [0,0], '--k')
legend('control', 'drink', 'inject')
xlabel('valve #')
ylabel('zscore')
hold off

%% figure 6
control = spkZ(contains(allClusters.group,'control'), :, :);
drink  = spkZ(contains(allClusters.group,'drink'), :, :);
inject = spkZ(contains(allClusters.group,'injected'), :, :);

figure(6)
clf

subplot(1,3,1)
hold on
valveGroup = 1;

y = squeeze(mean(control(:,:,valveGroup),1));
err = std(control(:,:,valveGroup),0,1) ./ sqrt(size(control,1));
shadedErrorBar(binEdges(1:end-1), y,err, 'lineProps', 'b')

y = squeeze(mean(drink(:,:,valveGroup),1));
err = std(drink(:,:,valveGroup),0,1) ./ sqrt(size(drink,1));
shadedErrorBar(binEdges(1:end-1), y,err, 'lineProps', 'r')

y = squeeze(mean(inject(:,:,valveGroup),1));
err = std(inject(:,:,valveGroup),0,1) ./ sqrt(size(inject,1));
shadedErrorBar(binEdges(1:end-1), y,err, 'lineProps', 'g')

ylabel('zscore')
title('valves 1-10')
xlabel('time (s)')
ylim([-1,1])

subplot(1,3,2)
hold on
valveGroup = 3;

y = squeeze(mean(control(:,:,valveGroup),1));
err = std(control(:,:,valveGroup),0,1) ./ sqrt(size(control,1));
shadedErrorBar(binEdges(1:end-1), y,err, 'lineProps', 'b')

y = squeeze(mean(drink(:,:,valveGroup),1));
err = std(drink(:,:,valveGroup),0,1) ./ sqrt(size(drink,1));
shadedErrorBar(binEdges(1:end-1), y,err, 'lineProps', 'r')

y = squeeze(mean(inject(:,:,valveGroup),1));
err = std(inject(:,:,valveGroup),0,1) ./ sqrt(size(inject,1));
shadedErrorBar(binEdges(1:end-1), y,err, 'lineProps', 'g')

title('valves 21-30')
xlabel('time (s)')
ylim([-1,1])

subplot(1,3,3)
hold on
valveGroup = 6;

y = squeeze(mean(control(:,:,valveGroup),1));
err = std(control(:,:,valveGroup),0,1) ./ sqrt(size(control,1));
shadedErrorBar(binEdges(1:end-1), y,err, 'lineProps', 'b')

y = squeeze(mean(drink(:,:,valveGroup),1));
err = std(drink(:,:,valveGroup),0,1) ./ sqrt(size(drink,1));
shadedErrorBar(binEdges(1:end-1), y,err, 'lineProps', 'r')

y = squeeze(mean(inject(:,:,valveGroup),1));
err = std(inject(:,:,valveGroup),0,1) ./ sqrt(size(inject,1));
shadedErrorBar(binEdges(1:end-1), y,err, 'lineProps', 'g')

title('valves 51-60')
xlabel('time (s)')
ylim([-1,1])

legend('control','drink','inject')






