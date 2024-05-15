%% Simulation parameters
tstep = 1; % Size of timestep in seconds (smaller gives more accurate simulation results but will take longer to run)
tMax = 60*12; % Total length of simulation in seconds (Ephys recording was 12 minutes)
tStop = 60*2; % Time at which to stop the EtOH injection in seconds (2 minute injection during ephys recording) 
d = 1200; % um2 / s diffusion constant in brain (Gonzalez, June 1998)
ek = 0.8/60; %Elimination rate constant (proportion eliminated / second) (lower bound estimate given blood flow; Gonzalez, June 1998)
ek = 2/60; %Elimination rate constant (proportion eliminated / second) (upper bound estimate given model; Gonzalez, June 1998)

%% Run the simulation 
[conc, t, dist] = runSim(tstep,tStop, tMax, d, ek); % Run the simulation 
realConc =  conc*1422;% Convert the normalized concentration to real units (1.8% v/v Ethanol = 1422 mg/dL)

%% Save data for later analysis
diffusion = struct;
diffusion.coeff = d;
diffusion.ek = ek;
diffusion.conc = squeeze(realConc(:,200,200:end))';
diffusion.dist = dist(200,200:end)';
diffusion.time = [0,t];
diffusion.conc = [zeros(size(diffusion.conc,1),1) , diffusion.conc]; %Add 0 concentration to time 0
save('diffusion.mat', "diffusion")

%% plot [EtOH] for different distances from the injection site
figure(2); clf
hold on
thisDist = 38; %385 um corresponds to the closest 10% of units to injection
plot(t, squeeze(realConc(:,200+thisDist,200)))

thisDist = 56; %560 um corresponds to median unit distance from injection
plot(t, squeeze(realConc(:,200+thisDist,200)))

thisDist = 77;%775 um corresponds to the furthest 10% of units to injection
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

%% plot peak [EtOH] vs. distance from injection
figure(4);clf

plot(dist(200,200:end),peakConc(200,200:end), 'LineWidth',2)
% xlim([0,2000])
xlim([200,900])
% y_old = ylim;
% ylim([5,y_old(2)])
% set(gca, 'YScale', 'log')
ylim([0,500])
xlabel('Distance from injection (um)')
ylabel('Peak ethanol concentration (mg/dL)')

%% movie of concentration over time
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

%% Simulation function
function [ConcentrationTimecourse, t, dist2D] = runSim(tstep, tStop, tMax, d, ek)
    mult=10; % Voxel size in microns. Smaller is more accurate but used a lot more RAM
    sigma = sqrt(2*d*tstep) / mult; % Standard deviation of Gaussian (related to the fundamental solution to Fick's 2nd law of diffusion)
    [xD,yD,zD] = meshgrid(1:401);  % Make a grid of indices for the brain cube
    dists = (sqrt((xD-200).^2 + (yD-200).^2 + (zD-200).^2)) * mult; % calculate the distance of each voxel from the center of the brain cube.
    concentration = zeros(401, 401, 401); % The 3D brain cube that tracks ethanol concentration

    t = tstep:tstep:tMax; % Array of times after injection
    nInj = sum(t<tStop); % number of timesteps occuring during the injection

    fullSphere = 5e8; % volume of entire injection = 0.5 uL
    injectionSphere = fullSphere/nInj; % volume injected per timestep
    r = (3*injectionSphere/4/pi) .^ (1/3); % radius of injection sphere
    isInjectedSphere = dists < r(1); % Logical values indicating which brain cube voxels occur inside the injection sphere

    ConcentrationTimecourse = nan(length(t),401,401); % Matrix to save concentration values for each timestep (only 2D to save RAM)
    figure(1); clf
    for i=1:length(t) % For each timestep
        if i<=nInj % If an injection should be occurping at this timestep
            concentration(isInjectedSphere) = concentration(isInjectedSphere) + 1; % Add ethanol to the injection sphere (normalized to 1)
        end
        concentration = imgaussfilt3(concentration, sigma, "Padding","replicate"); % Apply 3D gaussian filter to simulate diffusion of ethanol (Fick's 2nd law)
        concentration = concentration * exp(-ek*tstep); % Calculate [EtOH] after elimination
        imagesc(squeeze(concentration(:,:,200)));colorbar;% Display the concentration in the brain cube (2D slice)
        title([num2str(i),'/',num2str(length(t))]) % Display the current progress of the simulation
        ConcentrationTimecourse(i,:,:) = concentration(:,:,200); % Save a 2D snapshot of the [EtOH] during this timestep
        pause(0.01) % Briefly pause to give the the figure enough time to update
    end

    dist2D = squeeze(dists(200,:,:)); % Save a 2D snapshot of distance from brain cube center
end
