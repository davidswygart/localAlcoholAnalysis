function calcIsiViolations(clusters)
clusters.isiViolations = cellfun(@hillFormula, clusters.spikeTimes);
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

function hill = hillFormula(spk)
ISIt = 0.0015; % biological threshold for ISI violation
ISImin = 0; % minimum ISI threshold enforced by the data recording system used.
ISIs = diff(spk); % array of ISI
Tr = max(spk); % duration of recording
Ns = length(spk); % number of spikes
hill = sum(ISIs<ISIt)*Tr/2/(Ns^2)/(ISIt-ISImin);
end