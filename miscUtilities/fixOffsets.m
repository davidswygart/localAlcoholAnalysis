
matFiles = ls('*.mat');
nFiles = size(matFiles,1);

saveFolder = 'withOffset';
mkdir(saveFolder)

for i=1:nFiles
    load(matFiles(i,:));

    events.timestamp = events.timestamp - offset(i);

    subplot(2,1,1)
    plot(events.timestamp(events.line==8), events.state(events.line==8))
    xlim([0,3000])
    
    subplot(2,1,2)
    allspk = sort(cell2mat(clusterInfo.spikeTimes));
    plot(allspk,1:length(allspk))
    xlim([0,3000])

    save([saveFolder,filesep,matFiles(i,:)], "channel","clusterInfo","events","paramsFile","waveformsTime")
end


