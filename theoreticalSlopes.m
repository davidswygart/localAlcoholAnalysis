t = 1:60;


control = t*0;
drink = -1 ./ (1+ exp(3-t*.1)) + .05;
inject = fliplr(drink);

figure(99)
clf
plot(t,control, 'b')
hold on
plot(t,drink, 'r')
plot(t,inject, 'g')
legend('Control', 'Drink', 'Inject')

ylabel('Zscore')
xlabel('Time')
ylim([-1, .1])

