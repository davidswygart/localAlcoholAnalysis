tstep = 1;
tMax = 60*12;
tStop = 60*2;
d = 1200; % um2 / s diffusion constant in brain (Gonzalez, June 1998)
ek = 0.8/60; %Elimination rate constant (proportion eliminated / second)
ek = 2/60; %Elimination rate constant (proportion eliminated / second)
[conc, t, dist] = runSim(tstep,tStop, tMax, d, ek);
realConc =  conc*1422;% mg/dL (1.8% v/v)

%% Save data for later analysis
diffusion = struct;
diffusion.coeff = d;
diffusion.ek = ek;
diffusion.conc = squeeze(realConc(:,200,200:end))';
diffusion.dist = dist(200,200:end)';
diffusion.time = [0,t];
diffusion.conc = [zeros(size(diffusion.conc,1),1) , diffusion.conc]; %Add 0 concentration to time 0
save('diffusion.mat', "diffusion")
%% plot at different distances
figure(2); clf
hold on
thisDist = 38;%385 um
plot(t, squeeze(realConc(:,200+thisDist,200)))
thisDist = 56;%560 um
plot(t, squeeze(realConc(:,200+thisDist,200)))
thisDist = 77;%775 um
plot(t, squeeze(realConc(:,200+thisDist,200)))
plot([120,120],ylim,'r--')
legend('380 um', '560 um', '775 um')
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
histogram(goodClusters.distFromInj(contains(goodClusters.group,'inject')),12)
ylabel('Cluster count (inject group)')


hold on

yyaxis right
plot(dist(200,200:end),peakConc(200,200:end), 'LineWidth',2)
% xlim([0,2000])
xlim([200,900])
% y_old = ylim;
% ylim([5,y_old(2)])
% set(gca, 'YScale', 'log')
ylim([0,500])
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
histogram(goodClusters.distFromInj(contains(goodClusters.group,'inject')),12)
ylabel('Cluster count (inject group)')


hold on

yyaxis right
plot(dist(200,200:end),meanConc(200,200:end), 'LineWidth',2)
xlim([200,900])
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
    c.Label.String = "Concentration (mg/dL)";
    clim([1,1400])
    set(gca,'ColorScale','log')
    title([num2str(t(i)),'s'])
    pause(0.01)
end

%%
function [conc, t, dist2D] = runSim(tstep, tStop, tMax, d, ek)
    mult=10;
    sigma = sqrt(2*d*tstep) / mult;
    [xD,yD,zD] = meshgrid(1:401); 
    dists = (sqrt((xD-200).^2 + (yD-200).^2 + (zD-200).^2)) * mult;
    img = zeros(401, 401, 401);

    t = tstep:tstep:tMax;
    nInj = sum(t<tStop);

    fullSphere = 5e8; % volume of entire injection = 0.5 uL
    injectionSphere = fullSphere/nInj; % volume injected per timestep
    r = (3*injectionSphere/4/pi) .^ (1/3); % radius of injection sphere
    isInjectedSphere = dists < r(1);

    conc = nan(length(t),401,401);
    figure(1); clf
    for i=1:length(t)
        if i<=nInj
            img(isInjectedSphere) = img(isInjectedSphere) + 1;
        end
        img = imgaussfilt3(img, sigma, "Padding","replicate");
        img = img * exp(-ek*tstep);
        imagesc(squeeze(img(:,:,200)));colorbar;%clim([0,2])
        title([num2str(i),'/',num2str(length(t))])
        conc(i,:,:) = img(:,:,200);
        pause(0.01)
    end

    dist2D = squeeze(dists(200,:,:));
end
