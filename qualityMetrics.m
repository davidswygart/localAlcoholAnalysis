%% collect all spike times
[goodSpikes, muaSpikes, noiseSpikes] = loadSpikes();
goodSpikes = removeBelowRate(goodSpikes);
muaSpikes = removeBelowRate(muaSpikes);
noiseSpikes = removeBelowRate(noiseSpikes);


%% look at ISI
figure(1)
subplot(3,1,1)
hillGood = plotISI(goodSpikes);
title("good")
subplot(3,1,2)
hillMua = plotISI(muaSpikes);
title("MUA")
subplot(3,1,3)
hillNoise = plotISI(noiseSpikes);
title("noise")
xlabel('log10 ISI violations (Hill formula)')


%% look at presence
figure(2)
subplot(3,1,1)
presenceGood = plotPresence(goodSpikes);
title("good")
subplot(3,1,2)
presenceMua = plotPresence(muaSpikes);
title("MUA")
subplot(3,1,3)
presenceNoise = plotPresence(noiseSpikes);
title("noise")

%% threshold based on the above critera
keepGood = hillGood < .5 & presenceGood > .9;
keepMua = hillMua < .5 & presenceMua > .9;
keepNoise = hillNoise < .5 & presenceNoise > .9;
%print mean(keepGood)
fprintf('proportion of good: %3.2f \n', mean(keepGood))
fprintf('number of good: %4.0f \n\n', sum(keepGood))
fprintf('proportion of mua: %3.2f \n', mean(keepMua))
fprintf('number of mua: %4.0f \n\n', sum(keepMua))
fprintf('proportion of noise: %3.2f \n', mean(keepNoise))
fprintf('number of noise: %4.0f \n\n', sum(keepNoise))


%%
function [goodSpikes, muaSpikes, noiseSpikes] = loadSpikes()
matFiles = ls('*.mat');
nFiles = size(matFiles,1);
goodSpikes = [];
muaSpikes = [];
noiseSpikes = [];
    for i=1:nFiles
        %concatonate all spikes
        load(matFiles(i,:))
        goodSpikes = [goodSpikes; spkGood];
        muaSpikes = [muaSpikes; spkMUA];
        noiseSpikes = [noiseSpikes; spkNoise];
    end
end

function hill = plotISI(spk)
hill = cellfun(@calcISIHill, spk);
histogram(log10(hill+0.00001), 20)
ylabel('counts')
xlim([-5.2,5])
hold on
L = ylim();
plot([log10(0.5),log10(0.5)],[0,L(2)], '--r' )
hold off
end

function hill = calcISIHill(spk)
isi = diff(spk);
thresh = 0.0015;
violations = sum(isi < thresh);
t = max(spk);
n = length(spk);
hill = violations*t / (2*n^2*thresh);
end


function spk = removeBelowRate(spk)
haveEnoughSpikes = cellfun(@(x) mean(diff(x)), spk)<10;
spk = spk(haveEnoughSpikes);
end

function [spkCounts, binEdges] = binSpikes(spk)
binSize = 60; %bin width in seconds
maxTime = 47*60; %max time in seconds
binEdges = 0:binSize:maxTime;
nbins = length(binEdges)-1;
nclusters = length(spk);
spkCounts = nan(nclusters,nbins);
for u =1:length(spk)
    spkCounts(u,:) = histcounts(spk{u}, binEdges);
end
end

function p = plotPresence(spk)
[spkCounts, ~] = binSpikes(spk);
p = mean(spkCounts>0,2);
histogram(p,20)
ylabel('counts')
hold on
L = ylim();
plot([.9,.9],[0,L(2)], '--r' )
hold off
end
