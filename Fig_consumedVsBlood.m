g = load("C:\Users\david\OneDrive - Indiana University\localAlcohol\ephysData\group.mat");
group = g.group;
group.mlPerkg = 1000 * group.fluidConsumed ./ group.mouseWeight;

control = group(contains(group.group,'control'),:);
injected = group(contains(group.group, 'injected'),:);
drink = group(contains(group.group, 'drink'),:);

dat = struct();
dat.control = control.mlPerkg;
dat.injected = injected.mlPerkg;
dat.drink = drink.mlPerkg;

figure(2)
clf
subplot(2,1,1)
violinplot(dat, {},'ShowMean', true, 'ShowMedian', false, 'ShowBox',true)
ylabel('volume consumed (ml/kg)')
oldYlim = ylim();
ylim([0,oldYlim(2)])

%%
density = 1.006;
concentrationAlcohol = 0.2;
x = drink.mlPerkg * concentrationAlcohol * density;
y = drink.brainAlcohol;

subplot(2,1,2)
scatter(x,y)
xlabel('alcohol consumed (g/kg)')
ylabel('brain alcohol (mg%)')

oldXlim = xlim();
xlim([0,oldXlim(2)])
