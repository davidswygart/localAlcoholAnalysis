
matFiles = ls('*.mat');
nFiles = size(matFiles,1);

saveFolder = 'withQualityMetrics';
mkdir(saveFolder)

for i=1:nFiles
    load(matFiles(i,:))
    clusterInfo = calcIsiViolations(clusterInfo);
    clusterInfo = calcPresenceRatio(clusterInfo);
    save([saveFolder,filesep,matFiles(i,:)], "channel","clusterInfo","events","paramsFile","waveformsTime")
end