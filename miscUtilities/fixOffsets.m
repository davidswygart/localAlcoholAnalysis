
matFiles = ls('*.mat');
nFiles = size(matFiles,1);

saveFolder = 'withOffset';
%mkdir(saveFolder)


bpodAfterSpikes = nan(nFiles,1);
tLessThan0AfterOffset = nan(nFiles,1);
hasNegative = nan(nFiles,1);

for i=1:nFiles
    load(matFiles(i,:));
    matFiles(i,:)

    %events.timestamp = events.timestamp - offset(i);
    figure(1)
    plot(events.timestamp(1:200))

    figure(2)
    subplot(2,1,1)
    plot(events.timestamp(events.line==8), events.state(events.line==8))
    xlim([-500,4000])
    
    subplot(2,1,2)
    allspk = sort(cell2mat(clusterInfo.spikeTimes));
    plot(allspk,1:length(allspk))
    xlim([-500,4000])



    hasNegative(i) = min(events.timestamp) < 0;
    bpodAfterSpikes(i) = max(events.timestamp(events.line==8)) > max(allspk);

    min(allspk)

    newT = events.timestamp - offset(i);
    tLessThan0AfterOffset(i) = min(newT(events.line==8)) < 0;

    %save([saveFolder,filesep,matFiles(i,:)], "channel","clusterInfo","events","paramsFile","waveformsTime")
end


