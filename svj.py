# vz stochastic volatility with jumps#
from numpy import log, sqrt, exp, arange, pi
from scipy import *
from pylab import *

# SVJ Put Price
def svj(S,K,v,t,r,dy,l,mj,sj,tv,kv,sv,rho, callput):
# where
#           S: stock price
#           K: Strike Price (level)
#           v: current variance
#           t: maturity
#           r: interest rate
#           dy: dividend yield
#
#   PARAMETRIC INPUTS
#         l: jump arrival
#         mj: mean jump
#         sj: s.d. of jump
#         tv: unconditonal variance
#         kv: speed of mean reversion
#         sv: vol. of vol.
#         rh: cor. b/w ret. and vol.
#         callput = 1 for call and = 0 for put


# PARAMETERS USED IN OTHER FUNCTIONS CALLED FROM THIS ONE. THEY NEED
# TO BE CHANGED IN THE FUNCTIONS THEMSELVES

#test case:
#svj(7.85,20,0.73**2,0.5,0.0576,0.,-0.47,0.96,1.18,0.98,0,0) =11.81


   ep= 0.00001;
   Ib0= 100;
   TOL=1e-4;
   dw=.005;


   w=arange(ep,Ib0,dw);

   x=log(S); # get log of stock;


   u1=.5;
   u2=-.5;
   jbar = (1 + mj)**(1j* w)* exp(sj**2/2.* w* (1j - w));
   a  = 1j* w* (r - dy - l* mj) + l * (1 + mj)* (jbar -1);
   A0 = u1 * 1j * w - .5* w**2;
   A1 =  w * 1j * sv* rho - (kv-rho* sv);
   A2 = .5* sv**2;
   k  = (-( rho*sv*1j*w - (kv-rho*sv)) + sqrt((rho*sv*1j*w - (kv-rho*sv))**2  - \
   sv**2*(2* u1* 1j* w - w**2)) )/sv**2;
   
   f11=exp(1j*w* x + a* t + tv*(k*t - log(A1 + A2*(1 + \
   exp((A1 + 2*A2*k)*t))*k)/A2) + tv* log(A1+ 2* A2* k)/A2 - \
    (((-1 + exp((A1 + 2*A2*k)*t))*k*(A1 + A2*k))/ \
     (A1 + A2*(1 + exp((A1 + 2*A2*k)*t))*k))*v);

   f11mod=((exp(-1j*w*log(K))*f11)/(1j*w)).real;
   PI11=.5+1/pi*sum(f11mod*dw);



   jbar = (1 + mj)**(1j* w)* exp(-sj**2/2* w* (1j + w));
   a  = 1j* w* (r - dy - l* mj) + l*(jbar -1);
   A0 = u2* 1j* w - .5* w**2;
   A1 =  w* 1j* sv* rho - kv;
   A2 = .5* sv**2;
   k  = (-( rho*sv*1j*w - kv) + sqrt((rho*sv*1j*w - kv)**2  - \
   sv**2* (2* u2* 1j* w - w**2)) )/sv**2;

   f22=exp(1j* w* x + a* t + tv* (k*t - log(A1 + A2*(1 + \
   exp((A1 + 2*A2*k)*t))*k)/A2) + tv* log(A1+ 2* A2* k)/A2 - \
    (((-1 + exp((A1 + 2*A2*k)*t))*k*(A1 + A2*k))/ \
     (A1 + A2*(1 + exp((A1 + 2*A2*k)*t))*k))*v);

   f22mod=((exp(-1j*w*log(K))*f22)/(1j*w)).real;

   PI22=.5+1./pi*sum(f22mod*dw);

   if callput==1:
         
      F=exp(x)*exp(-dy*t)*PI11 - K*exp(-r*t)*PI22;
      Delta=PI11;
   elif callput==0:
      F=-exp(x)*exp(-dy*t)*(1-PI11) + K*exp(-r*t)*(1-PI22);
      Delta=PI11-1;
   else:
      disp('callput not specified')

   return array([F,Delta])

