T = readtable('5 series.xls');
x = table2array(T(:, 'X'));
y = table2array(T(:, 'Y'));
x = x(~isnan(x));
y = y(~isnan(y));
Q = trapz(x,y);
gamma=0.999;
N = 687;
ni = 0;
t = 0;
while ni <= 1 + t.*t
   ni = ni+1;
   t = tinv(gamma, ni);
end
a = min(y);
b = max(y);
rng default;
r = a + (b-a).*rand(1, N);
e = model_energy(Q, x, y, r);
histogram(e, ceil(N/ni));
%histogram(sqrt(2*abs(e)*1.60218e-11/9.10938356e-31), ceil(N/ni));
%xlabel('Величина вектора скорости (м\с)')
xlabel('Энергия электронов (эВ)')
ylabel('Число значений');
set(gca,'FontSize',20)