
load('C:\Users\dis006\OneDrive - Indiana University\localAlcohol\ephysData\group.mat')
matFiles = string(ls('*.mat'));
nFiles = size(matFiles,1);

allClusters = table();

for i = 1:nFiles
    d = load(matFiles(i));
    g = group(contains(group.matName, erase(matFiles(i),'.mat')), :);
    
    if isempty(g)
        continue
    end

    g = repmat(g, size(d.clusterInfo,1),1);
    clusters = d.clusterInfo;
    clusters = renamevars(clusters,'group','spikeInterfaceLabel');
    clusters.bpod(:) = extractBpodTimes(d.events);    
    clusters = [clusters, g];
    allClusters =[allClusters;clusters];
end