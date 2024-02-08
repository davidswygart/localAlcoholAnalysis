function clusters = calcIsiViolations(clusters)
clusters.isiViolations = cellfun(@hillFormula, clusters.spikeTimes);
end

function hill = hillFormula(spk)
ISIt = 0.0015; % biological threshold for ISI violation
ISImin = 0; % minimum ISI threshold enforced by the data recording system used.
ISIs = diff(spk); % array of ISI
Tr = max(spk); % duration of recording
Ns = length(spk); % number of spikes
hill = sum(ISIs<ISIt)*Tr/2/(Ns^2)/(ISIt-ISImin);
end