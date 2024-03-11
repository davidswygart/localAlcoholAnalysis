%load("allClusters.mat")
allClusters = allClusters(allClusters.fr>0.1, :);
allClusters = allClusters(~contains(allClusters.phyLabel,'noise'), :);
allClusters = allClusters(allClusters.presenceRatio>0.9, :);
allClusters = allClusters(allClusters.isiViolations<1, :);
%%
isControl = contains(allClusters.group,'control');
isDrink =  contains(allClusters.group,'drink');
isInject = contains(allClusters.group,'inject');
%% choose bin width
binWidth=0.1;
binEdges = -1:binWidth:14;

%%
spkCounts  = binAroundValve(allClusters, binEdges,'smooth');

%% Group by N valves 
valveBinning = 10;
clumped = nan(size(spkCounts,1),size(spkCounts,2),size(spkCounts,3)/valveBinning);
for i=1:size(spkCounts,3)/valveBinning
    ind1 = (i-1)*valveBinning+1;
    ind2 = ind1 + valveBinning - 1;
    clumped(:,:,i) = mean(spkCounts(:,:,ind1:ind2), 3);
end

spkZ = zscore(clumped,0,2);

%% figure 6 plot spike train for various clumps of valves
figure(6)
clf
x = binEdges(1:end-1);
limx = [-1,6];

subplot(2,2,1)
hold on
valveGroup = 1;
addShadedLine(x,spkZ(isControl,:,valveGroup),'b','Control')
addShadedLine(x,spkZ(isDrink,:,valveGroup),'r','Drink')
addShadedLine(x,spkZ(isInject,:,valveGroup),'g','Inject')
plot(xlim(), [0,0], '--k')
ylabel('zscore')
title('valves 1-10')
xlabel('time (s)')
xlim(limx)
hold off

subplot(2,2,2)
hold on
valveGroup = 2;
addShadedLine(x,spkZ(isControl,:,valveGroup),'b','Control')
addShadedLine(x,spkZ(isDrink,:,valveGroup),'r','Drink')
addShadedLine(x,spkZ(isInject,:,valveGroup),'g','Inject')
plot(xlim(), [0,0], '--k')
title('valves 11-20')
xlabel('time (s)')
xlim(limx)
hold off

subplot(2,2,3)
hold on
valveGroup = 4;
addShadedLine(x,spkZ(isControl,:,valveGroup),'b','Control')
addShadedLine(x,spkZ(isDrink,:,valveGroup),'r','Drink')
addShadedLine(x,spkZ(isInject,:,valveGroup),'g','Inject')
plot(xlim(), [0,0], '--k')
title('valves 31-40')
xlabel('time (s)')
xlim(limx)
hold off

subplot(2,2,4)
hold on
valveGroup = 6;
addShadedLine(x,spkZ(isControl,:,valveGroup),'b','Control')
addShadedLine(x,spkZ(isDrink,:,valveGroup),'r','Drink')
addShadedLine(x,spkZ(isInject,:,valveGroup),'g','Inject')
plot(xlim(), [0,0], '--k')
title('valves 51-60')
xlabel('time (s)')
xlim(limx)
legend()
hold off

%% Pull out peak Inh and Exc Values
excWindow = val2Ind(binEdges, [0.3, 0.9]); %time where zscore is >0.5 for any trace in Fig 4
exc = squeeze(mean(spkZ(:,excWindow(1):excWindow(2),:), 2));
inhWindow = val2Ind(binEdges, [1.6,3.2]); %time where zscore is <-0.5 for any trace in Fig 4
inh = squeeze(mean(spkZ(:,inhWindow(1):inhWindow(2),:), 2));
%% figure 5
figure(5)
clf

subplot(2,1,1)
hold on
title('"Exc" peak (0.2-1.0s after valve)')
addShadedLine([],exc(isControl,:),'b','Control')
addShadedLine([],exc(isDrink,:),'r','Drink')
addShadedLine([],exc(isInject,:),'g','Inject')
plot(xlim(), [0,0], '--k')
ylabel('zscore')
hold off

subplot(2,1,2)
hold on
title('"Inh" peak (2-2.5s after valve)')
addShadedLine([],inh(isControl,:),'b','Control')
addShadedLine([],inh(isDrink,:),'r','Drink')
addShadedLine([],inh(isInject,:),'g','Inject')
plot(xlim(), [0,0], '--k')
legend('control', 'drink', 'inject')
xlabel('valve #')
ylabel('zscore')
hold off
%% save for stats in R
saveCsvForR_repeatedMeasuresMixedModel( ...
    {exc(isControl,:),exc(isDrink,:),exc(isInject,:)}, ...
    {"control","drink","inject"}, ...
    'ExcPeak.csv')
saveCsvForR_repeatedMeasuresMixedModel( ...
    {inh(isControl,:),inh(isDrink,:),inh(isInject,:)}, ...
    {"control","drink","inject"}, ...
    'InhPeak.csv')


%%
function vals = val2Ind(x, vals)
    for i=1:length(vals)
        [~,minInd]=min(abs(x-vals(i)));
        vals(i)=minInd;
    end
end



