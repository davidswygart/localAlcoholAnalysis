dataDir = 'C:\Users\dis006\OneDrive - Indiana University\localAlcohol\ephysData\';
binSize = 0.1; %bin width in seconds
maxTime = 55*60; %max time in seconds
xA = 0:binSize:maxTime;

matFiles = ls([dataDir, '*.mat']);
nFiles = size(matFiles,1);

for i=1:nFiles
    load([dataDir, matFiles(i,:)])

    spk = [spkGood; spkMUA; spkNoise];
    

    nGood = length(spkGood);
    nMUA = length(spkMUA);
    nNoise = length(spkNoise);

    clusterIdentity = zeros([nMUA+nGood+nNoise, 1]);
    clusterIdentity(1:nGood) = 2;
    clusterIdentity(nGood+1: nGood+nMUA) = 1;
    
    nclusters = length(spk);
    nbins = length(xA)-1;
    spkCounts = nan(nclusters,nbins);
    spkSmoothed = nan(nclusters,nbins);

    for u =1:length(spk)
        spkCounts(u,:) = histcounts(spk{u}, xA);
        meanISI = mean(diff(spk{u}));

        if meanISI < 11
            spkSmoothed(u,:) = imgaussfilt(spkCounts(u,:),meanISI/binSize);
        else
            spkSmoothed(u,:) = nan;
        end
    end


    spkSmoothedNorm = spkSmoothed ./ max(spkSmoothed, [],2);
    spkZscore = zscore(spkSmoothed, 0, 2);


    figure(1)
    imagesc(spkCounts)
    colorbar
    title("spike counts")

    figure(2)
    imagesc(spkSmoothed)
    colorbar
    title("spikes smoothed")

    figure(3)
    imagesc(spkSmoothedNorm)
    colorbar
    title("spike smoothed normalized")


    figure(4)
    imagesc(spkZscore)
    colorbar
    title("spike smoothed Zscore")

% waveforms
%     ind = 1;
%     [~,bestInd] = max(max(abs(waveformsGood(ind,:,:)),[],2));
%     plot(waveformsTime,squeeze(waveformsGood(ind,bestInd,:)));

end