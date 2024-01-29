
%%
figure(1)

subplot(2,1,1)
plot(eventsTimestamp)
xlabel('timestamp #')
ylabel('timestamp value (s)')
hold on
plot([485000,491000,491000,485000,485000], [-4,-4,50,50,-4],'r')
hold off

subplot(2,1,2)
plot(eventsTimestamp)
xlabel('timestamp #')
ylabel('timestamp value (s)')
hold on
plot([485000,491000,491000,485000,485000], [-4,-4,50,50,-4],'r')
xlim([484500,491500])
ylim([-5,55])
hold off

%% which time delta is most common for each line

times = eventsTimestamp(eventsLine==2 & eventsState==1);
deltaT = diff(times);
1/mode(round(deltaT,4))

%{
Line identities
1: 40 Hz camera
2: 110 Hz camera
3: 30 Hz cameara
4: 1 Hz sync signal
8: bpod TTL
%}

%% How many times and how long are they stuck on -1
isStuck = eventsTimestamp==-1;
StuckCamera = isStuck(eventsLine==2 & eventsState==1);

stuckStart = find(diff(StuckCamera) == 1);
stuckStop = find(diff(StuckCamera) == -1);

figure(2)
subplot(2,1,2)
timeStuck = (stuckStop+1-stuckStart) / 109.8901;
histogram(timeStuck,15)
xlabel('time stuck at -1 (s)')
ylabel('counts')

subplot(2,1,1)
timeBetweenDrops = diff(stuckStart) / 109.8901;
histogram(timeBetweenDrops)
ylabel('counts')
xlabel('time between timestamp resets (s)')

%& eventsLine==2 & eventsState==1);
%plot(isStuck)




%% Where do breaks occur?
resetPointInds = find(diff(eventsTimestamp)<-0.5);
resumeNormalInds = find(diff(eventsTimestamp)>1.5);
fixedTimeStamps = eventsTimestamp;

aOnTimes = eventsTimestamp(eventsLine==1 & eventsState==1);
aOnInds = find(eventsLine==1 & eventsState==1);

for rpi = 1:length(resetPointInds)
    lastGoodAOn = aOnTimes(aOnInds<resetPointInds(rpi));
    lastGoodAOn = lastGoodAOn(end);

    NeedFixed = aOnInds>resetPointInds(rpi) & aOnInds<resumeNormalInds(rpi);
    start = lastGoodAOn+aOnDiff;
    stop = lastGoodAOn+aOnDiff*sum(NeedFixed);
    aOnTimes(NeedFixed) = start:aOnDiff:stop;
    
    remaining = find(NeedFixed,1,'last')+1;
    aOnTimes(remaining:end) = aOnTimes(remaining:end) - aOnTimes(remaining);
    aOnTimes(remaining:end) = aOnTimes(remaining:end) + stop+aOnDiff;

    %add final value to the rest of the train
end


%% Get average camera framerates (time between on events and delay from on event to off event)
[aOnDiff, aOffWait] = medianDiff(eventsTimestamp,eventsLine,eventsState, 1);
[bOnDiff, bOffWait] = medianDiff(eventsTimestamp,eventsLine,eventsState, 2);
[cOnDiff, cOffWait] = medianDiff(eventsTimestamp,eventsLine,eventsState, 3);
fixedTimeStamps = eventsTimestamp;

% Simple reconstruct, just use last good value for each camera
aOnInds = eventsLine==1 & eventsState==1;
% aOnVals = eventsTimestamp(aOnInds);
% fixedTimeStamps(aOnInds) = repairWithDelta(aOnVals, aOnDiff);
% aOffInds = eventsLine==1 & eventsState==0;
% aOffVals = eventsTimestamp(aOffInds);
% aOnVals = fixedTimeStamps(aOnInds);
% fixedTimeStamps(aOffInds) = repairOffUsingON(aOnVals, aOffVals, aOffWait);

bOnInds = eventsLine==2 & eventsState==1;
bOnVals = eventsTimestamp(bOnInds);
fixedTimeStamps(bOnInds) = repairWithDelta(bOnVals, bOnDiff);
bOffInds = eventsLine==2 & eventsState==0;
bOffVals = eventsTimestamp(bOffInds);
bOnVals = fixedTimeStamps(bOnInds);
fixedTimeStamps(bOffInds) = repairOffUsingON(bOnVals, bOffVals, bOffWait);

% cOnInds = eventsLine==3 & eventsState==1;
% cOnVals = eventsTimestamp(cOnInds);
% fixedTimeStamps(cOnInds) = repairWithDelta(cOnVals, cOnDiff);
% cOffInds = eventsLine==3 & eventsState==0;
% cOffVals = eventsTimestamp(cOffInds);
% cOnVals = fixedTimeStamps(cOnInds);
% fixedTimeStamps(cOffInds) = repairOffUsingON(cOnVals, cOffVals, cOffWait);

%% repair 
goodInds =  eventsLine==2 ;
%badInds = eventsLine==4 | eventsLine==8;

% if(any(~(goodInds | badInds)))
%     warning("I'm missing some inds somewhere")
% end

interpStamps = interp1(find(goodInds),fixedTimeStamps(goodInds), 1:length(eventsLine));


% check 1 Hz signal
plot(diff(interpStamps(eventsLine==4 & eventsState==1)))
%% functions
function [onDiff, offWait] = medianDiff(eventsTimestamp,eventsLine,eventsState, line)
onDiff = diff(eventsTimestamp(eventsLine==line & eventsState==1));
onDiff = onDiff(onDiff>0 & onDiff<1);
plot(onDiff)
onDiff = mean(onDiff);

offTimes = eventsTimestamp(eventsLine==line & eventsState==0);
onTimes = eventsTimestamp(eventsLine==line & eventsState==1);

if length(offTimes)<length(onTimes)
    onTimes = onTimes(1:end-1);
end
offWait = offTimes-onTimes;
offWait = offWait(offWait>0 & offWait<1);
plot(offWait)
offWait = mean(offWait);
end

function timestamps = repairWithDelta(timestamps, delta)
plot(timestamps)
lastGoodInd = find(timestamps<0,1)-1;
lastGoodVal = timestamps(lastGoodInd);
numToFix = length(timestamps(lastGoodInd:end))-1;

endVal = lastGoodVal+numToFix*delta;
timestamps(lastGoodInd:end) = lastGoodVal:delta:endVal;
end

function OffT = repairOffUsingON(OnT, OffT, delta)
firstBadInd = find(OffT==-1,1);
firstOnInd = find(OnT > OffT(firstBadInd-1),1);
corrected = OnT(firstOnInd:end) + delta;
if length(corrected) > length(OffT(firstBadInd:end))
    corrected = corrected(1:end-1);
end

OffT(firstBadInd:end) = corrected;
end

