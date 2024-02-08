
matFiles = ls('*.mat');
nFiles = size(matFiles,1);

for i=1:nFiles
    load(matFiles(i,:))
    clusterInfo = calcIsiViolations(clusterInfo);
    clusterInfo = calcPresenceRatio(clusterInfo);
end