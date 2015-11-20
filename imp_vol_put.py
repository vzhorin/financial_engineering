from scipy import stats
from numpy import log, sqrt, exp
from pylab import *
import time

def imp_vol_put(P,S,X,rd,q,T):

#f= imp_vol_put(P,S,X,rd,q,T) computes the implied volatility
# for a put when the traded call price is P.
# P and X must be column vectors of equal size

   P = P;
   X = X;
   T = T*1.;

   NP=len(P);
   NX=len(X);
 
   if NP == NX:
            
      S = S*1.;
      r = rd*1.;
      q = q*1.;

      ff=[];

      for np in range(NP):

         sigma=.3;
         error=.000001;

         dv=error+1;
         tic = time.time();
         
         while abs(dv) > error:
            d1=(log(S/X[np])+(r-q+sigma**2/2.)*T)/(sigma*sqrt(T));
            
            d2=d1-sigma*sqrt(T);
            nd1=stats.norm.cdf(-d1);
            nd2=stats.norm.cdf(-d2);
            npd1=stats.norm.pdf(d1);

            PriceError=-S*exp(-q*T)*nd1+X[np]*exp(-r*T)*nd2-P[np];
            Vega=S*sqrt(T)*exp(-q*T)*npd1;
            if Vega==0:
               disp('No Volatility can be found')
               sigma=NaN;
               break
            
            dv=PriceError/Vega;
            sigma=sigma-dv;
            time2=time.time()-tic;
            if time2>60:
               disp('the routine did not converge within 60 seconds')
               sigma=NaN;
               break
         ff.append(sigma);
   else:
      disp('P and X are not of equal size')
      ff=NaN;
      
   f=ff;

   return f



