binSize = 0.1; %bin width in seconds
maxTime = 50*60; %max time in seconds
subRefractoryThreshold = 0.002;
xA = 0:binSize:maxTime;
nbins = length(xA)-1;
trialEventLn = 8; % digital IN line# on DAQ board getting trial signal from BPOD

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


    %% get bpod ttl signal
%{
Description of bpod event encoding
1: baseline (10 minutes)
	on 1s
	off 1s
2: microinjection (2 minutes)
	on 0.5s
	off 1s
3: Post injection (10 minutes)
	on 2s
	off 1s
4: Sipper time (15 minutes)
	on: 4s pretrial
 	off: 11.1s (0.1 s valve open ->  10 s consumption -> 1s suction)
5: Tail time (10 minutes)
	on: 3s
	off: 1s
%}
    bpodTimeStamp = eventsTimestamp(eventsLine == trialEventLn);
    bpodState = eventsState(eventsLine == trialEventLn);
    



%%


if any(bpodState(1:2:end) ~= 1) || any(bpodState(2:2:end) ~= 0)
    error("bpod TTL state does not perfectly alternate between true-false")
end

bpodON = bpodTimeStamp(1:2:end);
bpodOFF = bpodTimeStamp(2:2:end);
onLength = round(bpodOFF-bpodON, 1);

baselineStart = bpodON(find(onLength==1, 1))/ binSize;
microInjectionStart = bpodON(find(onLength==0.5, 1))/ binSize;
postInjectionStart = bpodON(find(onLength==2, 1))/ binSize;
SipperStart = bpodON(find(onLength==4, 1))/ binSize;
tailStart = bpodON(find(onLength==3, 1))/ binSize;

if (~issorted([baselineStart, microInjectionStart, postInjectionStart, SipperStart, tailStart]))
    error("start times of bpod stages are not in the correct order")
end

dropDeployed = bpodOFF(onLength==4)/ binSize;

if abs(length(dropDeployed) - 60) > 1 %check of off by more than 1
    warning(['Expected 60 sipper drops. Only found ', num2str(length(dropDeployed))])
end



    %%
    figure(2)
    subplot(2,1,1)
    imagesc(spkZscore, [-2 2])
    colorbar
    title("spike Zscore")
    hold on
    firstNoiseInd = find(clusterIdentity==0, 1);
    plot([0, size(spkCounts,2)],[firstNoiseInd,firstNoiseInd], 'r')
    
    plot([baselineStart,baselineStart], [-20,0], 'r')
    text(baselineStart,-15,' Baseline')
    plot([microInjectionStart,microInjectionStart], [-20,0], 'r')
    text(microInjectionStart,-15,' Inj')
    plot([postInjectionStart,postInjectionStart], [-20,0], 'r')
    text(postInjectionStart,-15,' Post-inject')
    plot([SipperStart,SipperStart], [-20,0], 'r')
    text(SipperStart,-15,' Sipper')
    plot([tailStart,tailStart], [-20,0], 'r')
    text(tailStart,-15,' Tail')
    scatter(dropDeployed, ones(length(dropDeployed),1)*-1, '.r')

    currentYlim = ylim();
    ylim([-21,currentYlim(2)])
    xlabel('time (s*10)')
    ylabel('cluster #')
    hold off
    
    subplot(2,1,2)
    subRefractoryCounts(subRefractoryCounts==0) = nan;
    imagesc(subRefractoryCounts, 'alphadata', 0+~isnan(subRefractoryCounts))
    colorbar
    title(['Subrefractory spikes (ISI < ' num2str(subRefractoryThreshold) 's)'])
    hold on
    plot([0, size(spkCounts,2)],[firstNoiseInd,firstNoiseInd], 'r')
    
    plot([baselineStart,baselineStart], [-20,0], 'r')
    text(baselineStart,-15,' Baseline')
    plot([microInjectionStart,microInjectionStart], [-20,0], 'r')
    text(microInjectionStart,-15,' Inj')
    plot([postInjectionStart,postInjectionStart], [-20,0], 'r')
    text(postInjectionStart,-15,' Post-inject')
    plot([SipperStart,SipperStart], [-20,0], 'r')
    text(SipperStart,-15,' Sipper')
    plot([tailStart,tailStart], [-20,0], 'r')
    text(tailStart,-15,' Tail')
    scatter(dropDeployed, ones(length(dropDeployed),1)*-1, '.r')

    currentYlim = ylim();
    ylim([-21,currentYlim(2)])
    ylabel('cluster #')
    hold off
    

    %% percentage of subrefractory spikes
    proportionSubRefractory = sum(subRefractoryCounts,2, 'omitnan') ./ sum(spkCounts,2);

    figure(3)
    plot([firstNoiseInd,firstNoiseInd], [0,1], 'r')
    hold on
    text(firstNoiseInd,.7,' "Noise" clusters')
    bar(proportionSubRefractory)
    xlabel('cluster #')
    ylabel(['proportion of spikes with ISI < ' num2str(subRefractoryThreshold) 's'])
    hold off




% waveforms
%     ind = 1;
%     [~,bestInd] = max(max(abs(waveformsGood(ind,:,:)),[],2));
%     plot(waveformsTime,squeeze(waveformsGood(ind,bestInd,:)));

end


