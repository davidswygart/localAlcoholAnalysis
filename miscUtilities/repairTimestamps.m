allVars = who;

%%
figure(1)

subplot(2,1,1)
plot(events.timestamp)
xlabel('timestamp #')
ylabel('timestamp value (s)')
hold on
plot([485000,491000,491000,485000,485000], [-4,-4,50,50,-4],'r')
hold off

subplot(2,1,2)
plot(events.timestamp)
xlabel('timestamp #')
ylabel('timestamp value (s)')
hold on
plot([485000,491000,491000,485000,485000], [-4,-4,50,50,-4],'r')
xlim([484500,491500])
ylim([-5,55])
hold off

%% which time delta is most common for each line

times = events.timestamp(events.line==2 & events.state==1);
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
isStuck = events.timestamp==-1;
StuckCamera = isStuck(events.line==2 & events.state==1);

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

%& events.line==2 & events.state==1);
%plot(isStuck)




%% Where do breaks occur?
% resetPointInds = find(diff(events.timestamp)<-0.5);
% resumeNormalInds = find(diff(events.timestamp)>1.5);
% fixedTimeStamps = events.timestamp;
% 
% aOnTimes = events.timestamp(events.line==1 & events.state==1);
% aOnInds = find(events.line==1 & events.state==1);
% 
% for rpi = 1:length(resetPointInds)
%     lastGoodAOn = aOnTimes(aOnInds<resetPointInds(rpi));
%     lastGoodAOn = lastGoodAOn(end);
% 
%     NeedFixed = aOnInds>resetPointInds(rpi) & aOnInds<resumeNormalInds(rpi);
%     start = lastGoodAOn+aOnDiff;
%     stop = lastGoodAOn+aOnDiff*sum(NeedFixed);
%     aOnTimes(NeedFixed) = start:aOnDiff:stop;
% 
%     remaining = find(NeedFixed,1,'last')+1;
%     aOnTimes(remaining:end) = aOnTimes(remaining:end) - aOnTimes(remaining);
%     aOnTimes(remaining:end) = aOnTimes(remaining:end) + stop+aOnDiff;
% 
%     %add final value to the rest of the train
% end


%% Get average camera framerates (time between on events and delay from on event to off event)
[aOnDiff, aOffWait] = medianDiff(events, 1);
[bOnDiff, bOffWait] = medianDiff(events, 2);
[cOnDiff, cOffWait] = medianDiff(events, 3);
fixedTimeStamps = events.timestamp;

% Simple reconstruct, just use last good value for each camera
% aOnInds = events.line==1 & events.state==1;
% aOnVals = events.timestamp(aOnInds);
% fixedTimeStamps(aOnInds) = repairWithDelta(aOnVals, aOnDiff);
% aOffInds = events.line==1 & events.state==0;
% aOffVals = events.timestamp(aOffInds);
% aOnVals = fixedTimeStamps(aOnInds);
% fixedTimeStamps(aOffInds) = repairOffUsingON(aOnVals, aOffVals, aOffWait);

bOnInds = events.line==2 & events.state==1;
bOnVals = events.timestamp(bOnInds);
fixedTimeStamps(bOnInds) = repairWithDelta(bOnVals, bOnDiff);
bOffInds = events.line==2 & events.state==0;
bOffVals = events.timestamp(bOffInds);
bOnVals = fixedTimeStamps(bOnInds);
fixedTimeStamps(bOffInds) = repairOffUsingON(bOnVals, bOffVals, bOffWait);

% cOnInds = events.line==3 & events.state==1;
% cOnVals = events.timestamp(cOnInds);
% fixedTimeStamps(cOnInds) = repairWithDelta(cOnVals, cOnDiff);
% cOffInds = events.line==3 & events.state==0;
% cOffVals = events.timestamp(cOffInds);
% cOnVals = fixedTimeStamps(cOnInds);
% fixedTimeStamps(cOffInds) = repairOffUsingON(cOnVals, cOffVals, cOffWait);

%% repair 
goodInds =  events.line==2 ;
%badInds = events.line==4 | events.line==8;

% if(any(~(goodInds | badInds)))
%     warning("I'm missing some inds somewhere")
% end

interpStamps = interp1(find(goodInds),fixedTimeStamps(goodInds), 1:length(events.line));


% check 1 Hz signal
clf
plot(diff(interpStamps(events.line==4 & events.state==1)))

%%
events.timestamp = interpStamps';
save('repaired.mat', allVars{:})

%% functions
function [onDiff, offWait] = medianDiff(events, line)
onDiff = diff(events.timestamp(events.line==line & events.state==1));
onDiff = onDiff(onDiff>0 & onDiff<1);
plot(onDiff)
onDiff = mean(onDiff);

offTimes = events.timestamp(events.line==line & events.state==0);
onTimes = events.timestamp(events.line==line & events.state==1);

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

