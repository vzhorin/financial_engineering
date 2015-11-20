# vz stochastic volatility#
from numpy import log, sqrt, exp, arange, pi
from scipy import *
from pylab import *

# SV Call/Put Price
def sv(S,K,v,t,r,dy,lambd, tv,kv,sv,rho, callput):

# where
#           S: stock price
#           K: Strike Price (level)
#           v: current variance, vol^2
#           t: maturity
#           r: interest rate
#           dy: dividend yield
#
#   PARAMETRIC INPUTS
#         lambd - market price of risk  
#         tv: theta unconditonal variance
#         kv: kappa speed of mean reversion
#         sv: ss vol. of vol.
#         rh: cor. b/w ret. and vol.
#         callput = 1 for call and = 0 for put


#test case:
#sv(7.85,20,0.73**2,0.5,0.0576,0.,-0.47,0.96,1.18,0.98,0,0) =11.81

   ep = 0.00001;
   Ib0 = 100;
   TOL = 1e-4;
   dphi = .005;

   phi = arange(ep,Ib0, dphi);

   u1 = .5;
   u2 = -.5;
   b1 = kv + lambd - rho * sv;
   b2 = kv + lambd;
   a = kv * tv;

   d1 = sqrt((rho*sv*phi*1j-b1)**2-sv**2*(2*u1*phi*1j-phi**2));
   d2=sqrt((rho*sv*phi*1j-b2)**2-sv**2*(2*u2*phi*1j-phi**2));
   g1=(b1-rho*sv*phi*1j+d1)/(b1-rho*sv*phi*1j-d1);
   g2=(b2-rho*sv*phi*1j+d2)/(b2-rho*sv*phi*1j-d2);

   D1=(b1-rho*sv*phi*1j+d1)/sv**2*((1-exp(d1*t))/(1-g1*exp(d1*t)));
   D2=(b2-rho*sv*phi*1j+d2)/sv**2*((1-exp(d2*t))/(1-g2*exp(d2*t)));

   C1=(r-dy)*phi*1j*t+a/sv**2*((b1-rho*sv*phi*1j+d1)*t-2*log((1-g1*exp(d1*t))/(1-g1)));
   C2=(r-dy)*phi*1j*t+a/sv**2*((b2-rho*sv*phi*1j+d2)*t-2*log((1-g2*exp(d2*t))/(1-g2)));

   f1=exp(C1+D1*v+1j*phi*log(S));
   f2=exp(C2+D2*v+1j*phi*log(S));

   PI1=.5+1./pi*sum(((exp(-1j*phi*log(K))*f1)/(1j*phi)).real*dphi);
   PI2=.5+1./pi*sum(((exp(-1j*phi*log(K))*f2)/(1j*phi)).real*dphi);

   if callput == 1:
      F=S*exp(-dy*t)*PI1-K*exp(-r*t)*PI2;
      return F
   elif callput == 0:
      F=-S*exp(-dy*t)*(1-PI1) + K*exp(-r*t)*(1-PI2);
      return F
   else:
     disp('callput not specified')
     
   return
   
