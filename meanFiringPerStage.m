clear
matFiles = ls('*.mat');
nFiles = size(matFiles,1);


allClusters = table;

for i=1:nFiles
    load(matFiles(i,:))
    bpodTimes = extractBpodTimes(events);
    goodClusters = clusterInfo.fr>.1;
    goodClusters = goodClusters & ~contains(clusterInfo.phyLabel,'noise');
    goodClusters = goodClusters & clusterInfo.isiViolations<1;
    goodClusters = goodClusters & clusterInfo.presenceRatio>0.9;
    clusterInfo = clusterInfo(goodClusters,:);

    mfileName = erase(matFiles(i,:), '.mat')

    for c=1:size(clusterInfo,1)
        spk = clusterInfo.spikeTimes{c};
        clusterInfo.matName{c} = mfileName;

        cluseterInfo.stages{c} = {'baseline','injection','postInject','sipper','tail'};
        avgSpk = nan(1,5);
        avgSpk(1) = sum(spk>bpodTimes.baselineStart       & spk<bpodTimes.baselineStart+        60*10)/(60*10);
        avgSpk(2) = sum(spk>bpodTimes.microInjectionStart & spk<bpodTimes.microInjectionStart+  60* 2)/(60* 2);
        avgSpk(3) = sum(spk>bpodTimes.postInjectionStart  & spk<bpodTimes.postInjectionStart+   60*10)/(60*10);
        avgSpk(4) = sum(spk>bpodTimes.sipperStart         & spk<bpodTimes.sipperStart+          60*15)/(60*15);
        avgSpk(5) = sum(spk>bpodTimes.tailStart           & spk<bpodTimes.tailStart+            60*10)/(60*10); 


        clusterInfo.avgSpk(c) = {avgSpk};
    end
    if size(clusterInfo,1) > 0
        allClusters = [allClusters;clusterInfo];
    end
end
%%

control = allClusters(ismember(allClusters.matName, group.matName(contains(group.group, 'control'))),:);
drink = allClusters(ismember(allClusters.matName, group.matName(contains(group.group, 'drink'))),:);
inject = allClusters(ismember(allClusters.matName, group.matName(contains(group.group, 'injected'))),:);

avgControl = mean(cell2mat(control.avgSpk),1);
stdControl = std(cell2mat(control.avgSpk), [], 1);
nControl = size(control,1);

avgDrink= mean(cell2mat(drink.avgSpk),1);
stdDrink = std(cell2mat(drink.avgSpk), [], 1);
nDrink = size(drink,1);

avgInject = mean(cell2mat(inject.avgSpk),1);
stdInject = std(cell2mat(inject.avgSpk), [], 1);
nInject = size(inject,1);


names = ["baseline","inject","post-inj","sipper","post-sip"];
figure(1)
clf
X = categorical({'baseline','inject','post-inj','sipper','post-sip'});
X = reordercats(X,{'baseline','inject','post-inj','sipper','post-sip'});
bar(X,[avgControl',avgDrink',avgInject'])
ylabel('mean spike rate')
legend('Control','Drink alcohol', 'Inject alcohol')



figure(2)
clf
subplot(3,1,1)
violinplot(cell2mat(control.avgSpk),names)
ylim([0,15])
ylabel('Control (spike rate)')
subplot(3,1,2)
violinplot(cell2mat(drink.avgSpk), names)
ylim([0,15])
ylabel('Drink alcohol (spike rate)')
subplot(3,1,3)
violinplot(cell2mat(inject.avgSpk), names)
ylim([0,15])
ylabel('Inject alcohol (spike rate)')


figure(3)
clf
subplot(3,1,1)
y = cell2mat(control.avgSpk);
y = y(:, 2:end) ./ y(:, 1);
y = log(y);
violinplot(y,names(2:end))
ylabel('Control: log(norm. spike rate)')
subplot(3,1,2)
y = cell2mat(drink.avgSpk);
y = y(:, 2:end) ./ y(:, 1);
y = log(y);
y(y<-5000) = -4 ;
violinplot(y,names(2:end))
ylabel('Drink alcohol: log(norm. spike rate)')
subplot(3,1,3)
y = cell2mat(inject.avgSpk);
y = y(:, 2:end) ./ y(:, 1);
y = log(y);
violinplot(y,names(2:end))
ylabel('Inject alcohol: log(norm. spike rate)')