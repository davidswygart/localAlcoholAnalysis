
%%
plot(eventsTimestamp(eventsLine==2 & eventsState==1))
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
%% Which lines keep resetting to -1?
unique(eventsLine(eventsTimestamp == -1))
% all of them

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

resetPoints = find(diff(eventsTimestamp)<0);

