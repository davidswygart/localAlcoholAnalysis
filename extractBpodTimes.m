function s = extractBpodTimes(events)
s = struct(); % output struct
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

bpodLine = 8; % digital IN line# on DAQ board getting trial signal from BPOD
events = events(events.line == bpodLine, :);

if any(events.state(1:2:end) ~= 1) || any(events.state(2:2:end) ~= 0)
    error("bpod TTL state does not perfectly alternate between true-false")
end

bpodON = events.timestamp(1:2:end);
bpodOFF = events.timestamp(2:2:end);
onLength = round(bpodOFF-bpodON, 1);

s.baselineStart = bpodON(find(onLength==1, 1));
s.microInjectionStart = bpodON(find(onLength==0.5, 1));
s.postInjectionStart = bpodON(find(onLength==2, 1));
s.sipperStart = bpodON(find(onLength==4, 1));
s.tailStart = bpodON(find(onLength==3, 1));

if (~issorted([s.baselineStart, s.microInjectionStart, s.postInjectionStart, s.sipperStart, s.tailStart]))
    error("start times of bpod stages are not in the correct order")
end

s.dropDeployed = bpodOFF(onLength==4);

if abs(length(s.dropDeployed) - 60) > 1 %check of off by more than 1
    warning(['Expected 60 sipper drops. Only found ', num2str(length(dropDeployed))])
end
end