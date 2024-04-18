function p = twoProportionZtest(count1, count2)
n1 = sum(count1);
n2 = sum(count2);

p1 = count1 / n1;
p2 = count2 / n2;

p = (count1+count2) / (n1+n2);

z = (p1-p2) ./ sqrt( p.*(1-p) * (1/n1 + 1/n2) );
p = 2 * (1 - normcdf(abs(z)));
end