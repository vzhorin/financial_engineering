hw = load('rate.txt');

index_nonzero = find(hw(:,3) > 0);
hw = hw(index_nonzero, 3);

t_bill = -365./91.*log(1-91*(hw/100./360.));

b = [ones(size(t_bill)-1,1) t_bill(1:end-1)]\t_bill(2:end);
err = t_bill(2:end) - b(1) - b(2)*t_bill(1:end-1);

dt = 1/252.; 
vas_e = b(1)/dt;
vas_g = (1-b(2))/dt;
vas_n = std(err)/sqrt(dt);
disp(vas_e); disp(vas_g);disp(vas_n);

hw = load('coupon.txt');

coupon_rate = hw(:,1)/100.; 
maturity = hw(:,3); 
bond_price = hw(:,4)/100.; 

num_bonds = size(bond_price,1);
disc_z = zeros(num_bonds,1);
term_structure = zeros(num_bonds,1);

disc_z(1) = bond_price(1)/(1.+coupon_rate(1)/2.);
term_structure(1) = -log(disc_z(1))/maturity(1);

for i=2:num_bonds
    temp = sum(disc_z(1:i-1)); 
    c = coupon_rate(i)/2.;
    disc_z(i) = (bond_price(i) - c*temp)/(1.+c);
    term_structure(i) = -log(disc_z(i))/maturity(i);
end

plot(disc_z, '.'); 
figure;
plot(term_structure, '.');

interest_rate = t_bill(end);
x = fminsearch(@(x) sum((vasicek(x, maturity,interest_rate, vas_e, vas_g, vas_n)-disc_z).^2), 0.);

risk_price = x;
disp(risk_price);
  
disc_z_vas = vasicek(risk_price, maturity,interest_rate, vas_e, vas_g, vas_n);
term_structure_vas = -log(disc_z_vas)./maturity;

figure;
plot(term_structure, '*'); hold on;plot(term_structure_vas); hold off;
