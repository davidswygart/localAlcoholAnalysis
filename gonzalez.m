tstep = 5;
tMax = 600;
tStop = tMax;
d = 1200; % um2 / s diffusion constant in brain (Gonzalez, June 1998)
ek = 2/60; %Elimination rate constant (proportion eliminated / second)
ek = 0.8/60; %Elimination rate constant (proportion eliminated / second)
[withElimination, t, dist, slice_withElimination] = runSim(tstep,tStop, tMax, d, ek);
%[noElimination, t, dist, slice_noElimination] = runSim(tstep,tStop, tMax, d, 0);
tStop = tMax-60*3;
%[withElimination_delay, t, dist, slice_withElimination_delay] = runSim(tstep, tStop, tMax, d, ek);

%% plot with and without elimination
figure(1); clf
hold on
thisDist = 100;
%plot(t, squeeze(noElimination(:,200+thisDist,200)))
plot(t, squeeze(withElimination(:,200+thisDist,200)))
legend('No elimination', 'With elimination')
xlabel('Time (s)')
ylabel('Concentration (norm)')
title('1mm from injection')

%% plot ethanol accross space at end of experiment
gonzalez_um = [0	250	500	750	1000	1250	1500	1750	2000];
gonzalez_measured = [1	0.7934782609	0.5652173913	0.402173913	0.2282608696	0.1086956522	0.07608695652	0.05434782609	0.03260869565];
gonzalez_modeled = [1	0.3768115942	0.1014492754	0.02898550725	0.007246376812	0	0	0	0];
gonzalez_postDecap3 = [1	0.3942307692	0.2307692308	0.1538461538	0.1057692308	0.04807692308	0.02884615385	0	0];

figure(3); clf; hold on
set(0, 'DefaultLineLineWidth', 2);
x = dist(200:end,200);


plot(gonzalez_um, gonzalez_measured, 'DisplayName', 'Gonzalez: measured')

plot(gonzalez_um, gonzalez_modeled, 'DisplayName', 'Gonzalez: model')

y = slice_withElimination(end,200:end);
y = y/ max(y);
plot(x, y, 'DisplayName', 'My model: elimination')

plot(gonzalez_um, gonzalez_postDecap3, 'DisplayName', 'Gonzalez: model+delay')

y = slice_withElimination_delay(end,200:end);
y = y/ max(y);
plot(x, y,  'DisplayName', 'My model: elimination+delay')

% y = slice_noElimination(end,200:end);
% y = y/ max(y);
% plot(x, y)
xlim([0,2000])
xlabel('distance')
legend('Location','best')

%% Compare ek
tstep = 5;
tMax = 600;
tStop = tMax-30;
d = 1200; % um2 / s diffusion constant in brain (Gonzalez, June 1998)

ek = 2/60; %Elimination rate constant (proportion eliminated / second)
[elim2, t, dist, slice_elim2] = runSim(tstep,tStop, tMax, d, ek);

ek = 1.4/60; %Elimination rate constant (proportion eliminated / second)
[elim1p4, t, dist, slice_elim1p4] = runSim(tstep,tStop, tMax, d, ek);

ek = 0.8/60; %Elimination rate constant (proportion eliminated / second)
[elim0p8, t, dist, slice_elim0p8] = runSim(tstep,tStop, tMax, d, ek);


%% plot to make sure we have reached steady state at 1mm
figure(1); clf
hold on
thisDist = 100;
plot(t, squeeze(elim2(:,200+thisDist,200)))
plot(t, squeeze(elim1p4(:,200+thisDist,200)))
plot(t, squeeze(elim0p8(:,200+thisDist,200)))
legend('Ek=2', 'Ek=1.4', 'Ek=0.8')
xlabel('Time (s)')
ylabel('Concentration (norm)')
title('1mm from injection')

%% plot ethanol accross space at end of experiment
figure(3); clf; hold on
set(0, 'DefaultLineLineWidth', 2);
x = dist(200:end,200);

plot(gonzalez_um, gonzalez_measured, '-k', 'DisplayName', 'Gonzalez: measured')
plot(gonzalez_um, gonzalez_modeled, '--k', 'DisplayName', 'Gonzalez: model')

y = slice_elim2(end,200:end);
y = y/ max(y);
plot(x, y, 'DisplayName', 'Ek = 2.0/min')

y = slice_elim1p4(end,200:end);
y = y/ max(y);
plot(x, y,  'DisplayName', 'Ek = 1.4/min')

y = slice_elim0p8(end,200:end);
y = y/ max(y);
plot(x, y,  'DisplayName', 'Ek = 0.8/min')


xlim([0,1400])
xlabel('Distance (um)')
ylabel('Concentration (norm to peak)')
legend('Location','best')
title('Post-Decapitation delay = 30s')
%% plot at different distances
figure(4); clf
hold on
thisDist = 75;
plot(t, squeeze(withElimination(:,200+thisDist,200)))
thisDist = 100;
plot(t, squeeze(withElimination(:,200+thisDist,200)))
thisDist = 125;
plot(t, squeeze(withElimination(:,200+thisDist,200)))
plot([120,120],ylim,'r--')
legend('750 um', '1000 um', '1250 um')
xlabel('Time (s)')
ylabel('Concentration (mg/dL)')
title(['Diffusion coefficient  = ',num2str(d,3),' um2/s'])



%%
function [conc, t, dist2D, sliceSum] = runSim(tstep, tStop, tMax, d, ek)
    mult=10;
    sigma = sqrt(2*d*tstep) / mult;
    [xD,yD,zD] = meshgrid(1:401); 
    dists = (sqrt((xD-200).^2 + (yD-200).^2 + (zD-200).^2)) * mult;
    img = zeros(401, 401, 401);

    t = tstep:tstep:tMax;
%     nInj = length(tstep:tstep:tStop);
%     fullSphere = 5e8; % 0.5 uL
%     rMax = (3*fullSphere/4/pi)^(1/3);
%     startingSphere = fullSphere/sum(t<tStop);
%     rMin =  (3*startingSphere/4/pi)^(1/3);
%     r = linspace(rMin,rMax, sum(t<tStop));
    r = 270/2 * ones(length(t));



    conc = nan(length(t),401,401);
    sliceSum = nan(length(t), size(img,1));
    figure(1); clf
    for i=1:length(t)
        if t(i)<=tStop
            %injVal = sum(dists<rMax,'all') / sum(dists < r(i), 'all') / nInj;
            %img(dists < r(i)) = img(dists < r(i)) + injVal;
            img(dists < r(i)) = 1;
        end
        img = imgaussfilt3(img, sigma, "Padding","replicate");
        img = img * exp(-ek*tstep);
        imagesc(squeeze(img(:,:,200)));colorbar;
        title([num2str(i),'/',num2str(length(t))])
        conc(i,:,:) = img(:,:,200);
        sliceSum(i,:) = squeeze(sum(img,[1,2]));
        pause(0.01)
    end

    dist2D = squeeze(dists(200,:,:));
end