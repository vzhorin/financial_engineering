from scipy import stats
from numpy import log, sqrt, exp
from pylab import *
import time

def imp_vol_call(C,S,X,rd,q,T): 
# computes the implied volatility
# for a call when the traded call price is C.
# C and X must be column vectors of equal size

   C=C;
   X=X;
   T=T*1.;

   NC=len(C);
   NX=len(X);

   if NC == NX:
            
      S=S*1.;
      r=rd*1.;
      q=q*1.;

      ff=[];

      for nc in range(NC):
            
         sigma=.3;
         error=.000001;

         dv=error+1.;
         tic = time.time();
         while abs(dv)>error:

            d1=(log(S/X[nc])+(r-q+sigma**2/2.)*T)/(sigma*sqrt(T));
            d2=d1-sigma*sqrt(T);
            nd1 = stats.norm.cdf(d1);
            nd2 = stats.norm.cdf(d2);
            npd1 = stats.norm.pdf(d1);

            PriceError=S*exp(-q*T)*nd1-X[nc]*exp(-r*T)*nd2-C[nc];
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
      disp('C and X are not of equal size')
      ff=NaN;

   f=ff;
   return f




