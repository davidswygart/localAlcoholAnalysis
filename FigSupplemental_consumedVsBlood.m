%load("C:\Users\david\OneDrive - Indiana University\localAlcohol\ephysData\group.mat");
group.mlPerkg = 1000 * group.fluidConsumed ./ group.mouseWeight;

control = group(contains(group.group,'control'),:);
injected = group(contains(group.group, 'injected'),:);
drink = group(contains(group.group, 'drink'),:);

dat = struct();
dat.control = control.mlPerkg;
dat.injected = injected.mlPerkg;
dat.drink = drink.mlPerkg;

figure(7)
clf
subplot(2,1,1)
violinplot(dat, {},'ShowMean', true, 'ShowMedian', false, 'ShowBox',true)
ylabel('volume consumed (ml/kg)')
oldYlim = ylim();
ylim([0,oldYlim(2)+2])

%%
density = 1.006;
concentrationAlcohol = 0.2;
x = drink.mlPerkg * concentrationAlcohol * density;
y = drink.brainAlcohol;

subplot(2,1,2)
scatter(x,y,'red','filled')
xlabel('alcohol consumed (g/kg)')
ylabel('brain alcohol (mg%)')

oldXlim = xlim();
xlim([0,oldXlim(2)]+.1)

hold on



%%

x = x(~isnan(y));
y = y(~isnan(y));

p = polyfit(x,y,1);

fity = polyval(p,x);

hold on
newX = [min(x),max(x)];

plot(newX,polyval(p,newX), 'r--')

R_squared = corr(polyval(p,x),y)^2;
text(.5,40,['R^2 = ',num2str(R_squared,2)])
hold off

saveas(gcf, 'consumedVsBlood.svg')

%%



%%
x = group.brainAlcohol;
x(isnan(x)) = 0;
y = group.numGoodClusters;
figure(666); clf
scatter(x,y,'k','filled')
xlabel('brain alcohol (mg%)')
ylabel('number of "good" clusters')
p = polyfit(x,y,1);
newX = [0,80];
fitY = polyval(p,newX);

hold on
plot(newX,fitY, 'r--')

R_squared = corr(polyval(p,x),y)^2;
text(40,30,['R^2 = ',num2str(R_squared,2)])
hold off


%%
% group.numGoodClusters(:) = nan;
% for r = 1:size(group,1)
%     group.numGoodClusters(r) = sum(contains(goodClusters.matName, group.matName(r)));
% end