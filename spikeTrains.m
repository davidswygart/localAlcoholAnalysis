binSize = 0.1; %bin width in seconds
maxTime = 50*60; %max time in seconds
subRefractoryThreshold = 0.002;
xA = 0:binSize:maxTime;
nbins = length(xA)-1;

matFiles = ls('*.mat');
nFiles = size(matFiles,1);

for i=1:nFiles
    %concatonate all spikes
    load(matFiles(i,:))

    spk = [spkGood; spkMUA; spkNoise];
    nGood = length(spkGood);
    nMUA = length(spkMUA);
    nNoise = length(spkNoise);
    clusterIdentity = zeros([nMUA+nGood+nNoise, 1]);
    clusterIdentity(1:nGood) = 2;
    clusterIdentity(nGood+1: nGood+nMUA) = 1;
    
    %Get rid of cluster with mean ISI > 10s
    haveEnoughSpikes = cellfun(@(x) mean(diff(x)), spk)<10;
    spk = spk(haveEnoughSpikes);
    clusterIdentity = clusterIdentity(haveEnoughSpikes);
    nclusters = length(spk);
    
    %Examine clusters with unusually low ISI
    figure(1)
    cID = 1;
    exampleSpk = spk{cID};
    plot(exampleSpk(2:end),log10(diff(exampleSpk)))
    xlabel('time (s)')
    ylabel('log10(ISI) (s)') 
    
    
    spkCounts = nan(nclusters,nbins);
    subRefractoryCounts = nan(nclusters,nbins);

    for u =1:length(spk)
        spkCounts(u,:) = histcounts(spk{u}, xA);
        
        
        isSubRefractory = diff(spk{u}) < subRefractoryThreshold;
        subSpikes = spk{u};
        subSpikes  = subSpikes(2:end);
        subSpikes = subSpikes(isSubRefractory);
        subRefractoryCounts(u,:) = histcounts(subSpikes, xA);
        
%         meanISI = mean(diff(spk{u}));
%         if meanISI < 11
%             spkSmoothed(u,:) = imgaussfilt(spkCounts(u,:),meanISI/binSize/4);
%         else
%             spkSmoothed(u,:) = nan;
%         end
    end

    
    


%     spkNorm = spkCounts ./ max(spkCounts, [],2);
%     figure(1)
%     imagesc(spkNorm)
%     colorbar
%     title("spike normalized")
%     hold on
%     plot([0, size(spkCounts,2)],[nGood+nMUA+1,nGood+nMUA+1], 'r')
%     hold off

    spkZscore = zscore(spkCounts, 0, 2);

    figure(2)
    subplot(2,1,1)
    imagesc(spkZscore, [-2 2])
    colorbar
    title("spike Zscore")
    hold on
    firstNoiseInd = find(clusterIdentity==0, 1);
    plot([0, size(spkCounts,2)],[firstNoiseInd,firstNoiseInd], 'r')
    hold off
    
    subplot(2,1,2)
    subRefractoryCounts(subRefractoryCounts==0) = nan;
    imagesc(subRefractoryCounts, 'alphadata', 0+~isnan(subRefractoryCounts))
    colorbar
    title('Subrefractory spikes')
    hold on
    plot([0, size(spkCounts,2)],[firstNoiseInd,firstNoiseInd], 'r')
    hold off
    

% waveforms
%     ind = 1;
%     [~,bestInd] = max(max(abs(waveformsGood(ind,:,:)),[],2));
%     plot(waveformsTime,squeeze(waveformsGood(ind,bestInd,:)));

end