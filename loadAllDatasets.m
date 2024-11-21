
paths = dataAndFigDirectoryPaths();
load([paths.data, 'group.mat'])
matFiles = dir([paths.individualRecordings,'*.mat']);
nFiles = length(matFiles);

allClusters = table();

for i = 1:nFiles
    name = matFiles(i).name;
    d = load([paths.individualRecordings, name]);
    isInTable = contains(group.matName, erase(name,'.mat'));
    g = group(isInTable, :);
    
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

allClusters = addQualityMetrics(allClusters);
allClusters = calcInjectionDist(allClusters);
%% choose good clusters
goodPhy = ~contains(allClusters.phyLabel,'noise');
goodRate = allClusters.fr > 0.1;
goodISI = allClusters.isiViolations < 1;
goodPR =  allClusters.presenceRatio > 0.9;
goodAC = allClusters.ampCutoff < 0.1; 
goodDepth = allClusters.depth<1100;
goodClusters = allClusters(goodPhy & goodRate & goodISI & goodPR & goodAC & goodDepth, :);

%% save clusters
save([paths.data,'allClusters.mat'],'allClusters')
save([paths.data,'goodClusters.mat'],'goodClusters')