# vz- basic BS put/call option pricing subroutines #
from scipy import stats
from numpy import log, sqrt, exp, array

def BSC (S, X, r, q, sigma, T):
	"""
	[C,Delta,Gamma,Theta,Vega,Rho]=BSC(S,X,r,q,sigma,T)
	Compute Black and Scholes Call option premium
	"""
	d1=(log(S/X)+(r-q+sigma**2/2)*T)/(sigma*sqrt(T));
	d2=d1-sigma*sqrt(T);

	nd1=stats.norm.cdf(d1);
	nd2=stats.norm.cdf(d2);
	npd1=stats.norm.pdf(d1);

	C=S*exp(-q*T)*nd1-X*exp(-r*T)*nd2;
	Delta=exp(-q*T)*nd1;
	Gamma=exp(-q*T)/(sigma*S*sqrt(T))*npd1;
	Theta=-sigma*S*exp(-q*T)*npd1/(2*sqrt(T))+q*S*nd1*exp(-q*T)-r*X*exp(-r*T)*nd2;
	Vega=S*sqrt(T)*exp(-q*T)*npd1;
	Rho=X*T*exp(-r*T)*nd2;
	return array([C, Delta, Gamma, Theta, Vega, Rho])
	
def BSP(S,X,r,q,sigma,T):
	"""
	[P,Delta,Gamma,Theta,Vega,Rho]=BSP(S,X,r,q,sigma,T)
	Compute Black and Scholes Put option premium
	"""
	d1=(log(S/X)+(r-q+sigma**2/2)*T)/(sigma*sqrt(T));
	d2=d1-sigma*sqrt(T);

	nd1=stats.norm.cdf(-d1);
	nd2=stats.norm.cdf(-d2);
	npd1=stats.norm.pdf(d1);

	P=-S*exp(-q*T)*nd1+X*exp(-r*T)*nd2;
	Delta=-exp(-q*T)*nd1;
	Gamma=exp(-q*T)/(sigma*S*sqrt(T))*npd1;
	Theta=-sigma*S*exp(-q*T)*npd1/(2*sqrt(T))-q*S*nd1*exp(-q*T)+r*X*exp(-r*T)*nd2;
	Vega=S*sqrt(T)*exp(-q*T)*npd1;
	Rho=-X*T*exp(-r*T)*nd2;
	return array([P,Delta,Gamma,Theta,Vega,Rho]);
