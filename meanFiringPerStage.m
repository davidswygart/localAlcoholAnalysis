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

%control.avgSpk