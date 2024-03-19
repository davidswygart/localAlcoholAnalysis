
bubble_cubicMicrons = 5e8; % 0.5 uL
rMin = 108/2; %% ID of 33G needle
rMax = (3*bubble_cubicMicrons/4/pi)^(1/3);

tstep = 5;
tMax = 60*12;
tStop = 60*2;
d = 1600; % um2 / s (free in water)
%d = 160;
[conc, t] = runSim(tstep,tStop, tMax, d, rMin, rMax);

%% plot

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
title('Diffusion coefficient  = ',num2str(d),' um2/s')
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
function [conc, t] = runSim(tstep, tStop, tMax, d, rMin, rMax)
    bubble_cubicMicrons = 5e8; % 0.5 uL
    rFull = (3*bubble_cubicMicrons/4/pi)^(1/3);

    mult=10;
    sigma = sqrt(2*d*tstep) / mult;
    [xD,yD,zD] = meshgrid(1:401); 
    dists = (sqrt((xD-200).^2 + (yD-200).^2 + (zD-200).^2)) * mult;
    img = zeros(401, 401, 401);

    t = tstep:tstep:tMax;
    nInj = length(tstep:tstep:tStop);

    r = linspace(rMin,rMax, length(t));
    conc = nan(length(t),401,401);

    figure(1)
    for i=1:length(t)
        if t(i)<tStop
            injVal = sum(dists<rFull,'all') / sum(dists < r(i), 'all') / nInj;
            img(dists < r(i)) = img(dists < r(i)) + injVal;
        end
        img = imgaussfilt3(img, sigma, "Padding","replicate");
        imagesc(squeeze(img(:,:,200)));colorbar;
        title([num2str(i),'/',num2str(length(t))])
        pause(0.01)
        conc(i,:,:) = img(:,:,200);
    end
end



