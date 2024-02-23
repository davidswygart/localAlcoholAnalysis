
load('C:\Users\david\OneDrive - Indiana University\localAlcohol\ephysData\group.mat')
matFiles = string(ls('*.mat'));
nFiles = size(matFiles,1);

allClusters = table();

for i = 1:nFiles
    d = load(matFiles(i));
    clusters = d.clusterInfo;
    clusters = renamevars(clusters,'group','spikeInterfaceLabel');

    clusters.bpod(:) = extractBpodTimes(d.events);

    g = group(contains(group.matName, erase(matFiles(i),'.mat')), :);
    g = repmat(g, size(d.clusterInfo,1),1);
    
    clusters = [clusters, g];
    allClusters =[allClusters;clusters];
end