tstep = 5;
tMax = 60*12;
tStop = 60*2;
d = 1600/3; % um2 / s (free in water)
[conc, t, dist] = runSim(tstep,tStop, tMax, d);

%% plot at different distances
realConc =  conc*1436;
figure(2); clf
hold on
thisDist = 75;
plot(t, squeeze(realConc(:,200+thisDist,200)))
thisDist = 100;
plot(t, squeeze(realConc(:,200+thisDist,200)))
thisDist = 125;
plot(t, squeeze(realConc(:,200+thisDist,200)))
plot([120,120],ylim,'r--')
legend('750 um', '1000 um', '1250 um')
xlabel('Time (s)')
ylabel('Concentration (mg/dL)')
title(['Diffusion coefficient  = ',num2str(d,3),' um2/s'])

%% image peak alcohol concentration
figure(3);clf
peakConc = squeeze(max(realConc,[],1));
imagesc(peakConc)
c = colorbar;
c.Label.String = 'Ethanol peak concentration (mg/dL)';

%% plot peak alcohol concentration
figure(4);clf

yyaxis left
histogram(goodClusters.distFromInj(contains(goodClusters.group,'inject')))
ylabel('Cluster count (inject group)')


hold on

yyaxis right
plot(dist(200,200:end),peakConc(200,200:end), 'LineWidth',2)
% xlim([0,2000])
xlim([550,1450])
% y_old = ylim;
% ylim([5,y_old(2)])
% set(gca, 'YScale', 'log')
ylim([0,300])
xlabel('Distance from injection (um)')
ylabel('Peak ethanol concentration (mg/dL)')

%% image avg alcohol concentration
figure(3);clf
meanConc = squeeze(mean(realConc,1));
imagesc(meanConc)
c = colorbar;
c.Label.String = 'Ethanol mean concentration (mg/dL)';

%% plot avg alcohol concentration
figure(4);clf

yyaxis left
histogram(goodClusters.distFromInj(contains(goodClusters.group,'inject')))
ylabel('Cluster count (inject group)')


hold on

yyaxis right
plot(dist(200,200:end),meanConc(200,200:end), 'LineWidth',2)
xlim([0,2000])
xlim([550,1450])
% y_old = ylim;
% ylim([1,y_old(2)])
%set(gca, 'YScale', 'log')
ylim([0,150])
xlabel('Distance from injection (um)')
ylabel('Average ethanol concentration (mg/dL)')
%% movie
figure(1);clf
for i=1:size(realConc,1)
    imagesc(squeeze(realConc(i,:,:)))
    xlabel('Distance (um / 10)')
    ylabel('Distance (um / 10)')
    c = colorbar;
    c.Label.String = "concentration (mg/dL)";
    clim([1,1436])
    set(gca,'ColorScale','log')
    title([num2str(t(i)),'s'])
    pause(0.02)
end

%%
function [conc, t, dist2D] = runSim(tstep, tStop, tMax, d)
    mult=10;
    sigma = sqrt(2*d*tstep) / mult;
    [xD,yD,zD] = meshgrid(1:401); 
    dists = (sqrt((xD-200).^2 + (yD-200).^2 + (zD-200).^2)) * mult;
    img = zeros(401, 401, 401);

    t = tstep:tstep:tMax;
    nInj = length(tstep:tstep:tStop);

    fullSphere = 5e8; % 0.5 uL
    rMax = (3*fullSphere/4/pi)^(1/3);
    startingSphere = fullSphere/sum(t<tStop);
    rMin =  (3*startingSphere/4/pi)^(1/3);
    r = linspace(rMin,rMax, sum(t<tStop));

    conc = nan(length(t),401,401);
    figure(1); clf
    for i=1:length(t)
        if t(i)<tStop
            injVal = sum(dists<rMax,'all') / sum(dists < r(i), 'all') / nInj;
            img(dists < r(i)) = img(dists < r(i)) + injVal;
        end
        img = imgaussfilt3(img, sigma, "Padding","replicate");
        imagesc(squeeze(img(:,:,200)));colorbar;
        title([num2str(i),'/',num2str(length(t))])
        pause(0.01)
        conc(i,:,:) = img(:,:,200);
    end

    dist2D = squeeze(dists(200,:,:));
end