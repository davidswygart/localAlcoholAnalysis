plot(diff(eventsTimestamp(eventsLine==1)))
%%
plot(eventsTimestamp(eventsLine==2 & eventsState==1))

%% which time delta is most common for each line

times = eventsTimestamp(eventsLine==6 & eventsState==1);
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
%% Which lines keep resetting to -1?
unique(eventsLine(eventsTimestamp == -1))
% all of them

%% How many times and how long are they stuck on -1
isStuck = eventsTimestamp==-1;
stuckStart = diff(isStuck) == 1;
stuckStop = diff(isStuck) == -1;
timesStuck = sum(stuckStart)

lengthStuck = nan(timeStuck,1);
for i=1:timeStuck
    stuckStart(i):stuckStart
end


%& eventsLine==2 & eventsState==1);
%plot(isStuck)