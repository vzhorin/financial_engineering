function discount_z = vasicek(risk_price, maturity, interest_rate, vas_e, vas_g, vas_n)

vas_e_adj = vas_e - risk_price*vas_n; 
%disp(risk_price);
%disp(vas_e_adj);

vas_B = (1. - exp(-vas_g*maturity))/vas_g;
vas_A = 1./vas_g^2*(vas_B - maturity)*(vas_e_adj*vas_g-vas_n^2/2.) - (vas_n^2*vas_B.^2)/(4.*vas_g);

discount_z = exp(vas_A - vas_B*interest_rate); 

