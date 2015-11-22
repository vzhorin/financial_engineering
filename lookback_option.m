%vz - convertable security
%price lookback option using a Monte Carlo simulation and the Longstaff-Schwartz method
%%
price_zero = 5;   
price_cap = 20;
price_cap_initial = 10;
nominal = 1000.;
coupon = .08;
%volatility
sig = .20;	
time_to_expiry = 5.;           
interest_rate = log(1+4.50/100/4.)*4.;
lookback_window = 16;
%can only lookback after 6 months
lookback_lock = 120;
seed = 2590834526890;

%time step - one trading day
dt = 1/252.;
%number of Monte-Carlo paths
num_path = 100;

time_period = time_to_expiry;
num_grid = floor(time_period/dt);

%conversion price discounts - from covenant
conv_disc = zeros(num_grid,1);    
for i = 7:15
    conv_disc(floor(i*21/dt/252:(i+1)*21/dt/252)) = 0.01*i;
end 
conv_disc((i+1)*21/dt/252:end) = 0.01*i;

%not discounts before 15 month  - from covenant
conv_dis(1:14*21/dt/252) =0;
conv_disc = 1- conv_disc; 
randn('state',seed);
delta_w = sqrt(dt)*randn(num_path,num_grid);
price_sim = price_zero*ones(num_path,num_grid);
price_adj = zeros(num_path,1);
%%
icount =1;
for i=2:num_grid
    delta_grid = delta_w(:,i);
	price_sim(:,i) = (1 + interest_rate*dt)*price_sim(:,i-1) + sig*price_sim(:, i-1).*delta_grid;
end

price_adj = min(price_sim(:,(end-lookback_window):end), [],2)*conv_disc(end);
number_shares = nominal./min(price_adj, price_cap);

principal = number_shares.*price_sim(:,end);
%discounted = exp(-interest_rate*time_to_expiry)*principal;

total = exp((coupon)*time_to_expiry)*principal;

payoff = NaN*ones(num_path,num_grid);
payoff(:,num_grid) = total;

%% starting MC simulation processing
for i = 1:num_grid-1
    grid_point = num_grid-i;
% only search where adjusted price with 22 days min is lower than current
% price discounted (to be implemented)
    price_adj = min(price_sim(:,max(1,grid_point-lookback_window):grid_point), [],2)*conv_disc(grid_point);
    if(grid_point <120)
       price_adj = price_cap_initial*ones(length(price_adj),1);
    end
    %in_the_money = find(price_adj-price_sim(:, grid_point)*conv_disc(grid_point) < 0); 
    %in_the_money = find(payoff(grid_point+1) > 0);
    %count_itm = length(in_the_money);
   
% compare payoffs from continuation and immediate exercise
   
	number_shares = nominal./min(price_adj, price_cap);

    principal = number_shares.*price_sim(:,grid_point);
    total = exp((coupon)*time_to_expiry*grid_point/num_grid)*principal;
    payoff(:,grid_point) = total;
   
 %discount all back   
    if( grid_point == num_grid-1) 
        y = (ones(num_path,1)*exp(-interest_rate*[1:num_grid-grid_point]*dt)).*payoff(:,grid_point+1:num_grid);
    else
        y = sum(((ones(num_path,1)*exp(-interest_rate*[1:num_grid-grid_point]*dt)).*payoff(:,grid_point+1:num_grid))')';
    end
    stock_price = price_sim(:,grid_point);
    
    %LS regression
    x = [ones(num_path,1), stock_price, stock_price.^2];
    [U,W,V] = svd(x);
    b = V*(W\(U'*y));

   plot(stock_price,x*b,'.',stock_price,payoff(:,grid_point),':')
   legend('Expected Payoff if Wait','Payoff Today if Exercise')
   xlabel('Stock Price')
   title('Estimation of Exercise Frontier')
   pause(.0001)
   
   stop_rule = find(payoff(:,grid_point) >= x*b);
   continue_rule = setdiff([1:num_path],stop_rule);

%payoff now, no payoff later
   payoff(stop_rule,grid_point+1:num_grid) = zeros(length(stop_rule),num_grid-grid_point);
%can get payoff only once
   payoff(continue_rule,grid_point) = zeros(length(continue_rule),1);
end
 
%%expected value
final = sum(((ones(num_path,1)*exp(-interest_rate*[2:num_grid]*dt)).*payoff(:,2:num_grid))')';

total = mean(final);
std_err = std(final)/sqrt(num_path);
 
disp(total); disp(std_err);