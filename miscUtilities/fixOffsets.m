
matFiles = ls('*.mat');
nFiles = size(matFiles,1);

saveFolder = 'withOffset';
%mkdir(saveFolder)


bpodTooEarly = [];
bpodTooLate = [];
eventsHasNegative = [];

for i=1:nFiles
    load(matFiles(i,:));
    

    %events.timestamp = events.timestamp - offset(i);
    figure(1)
    plot(events.timestamp(1:200))

    figure(2)
    subplot(2,1,1)
    plot(events.timestamp(events.line==8), events.state(events.line==8))
    xlim([-500,4000])
    
    subplot(2,1,2)
    allspk = sort(cell2mat(clusterInfo.spikeTimes));
    plot(allspk,1:length(allspk))
    xlim([-500,4000])

    if max(allspk) < max(events.timestamp(events.line==8))
        warning('bpod events extend beyond spike times')
        bpodTooLate = [bpodTooLate; matFiles(i,:)];
    end

    if min(allspk) > min(events.timestamp(events.line==8))
        warning('bpod events start before spikes')
        bpodTooEarly = [bpodTooEarly; matFiles(i,:)];
    end

    if min(events.timestamp) < 0
        warning(['min event timestamp = ', num2str(min(events.timestamp))])
        eventsHasNegative = [eventsHasNegative;matFiles(i,:)];
    end



    %save([saveFolder,filesep,matFiles(i,:)], "channel","clusterInfo","events","paramsFile","waveformsTime")
end


