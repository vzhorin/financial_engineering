function [PriceT, YieldT] = CIRFun(vec,r0,T,logFlag)

% [PriceT, YieldT] = CIRFun(vec,r0,T,h)
% Compute the price and yield curve of the CIR model:
%
%		 dr = (eta - gamma r) dt + sqrt(alpha r) dX
%
%   gamma=vec(1); eta=vec(2); alpha=vec(3); 
% B=2*(exp(phi1*T)-1)./((gamma+phi1)*(exp(phi1*T)-1)+2*phi1);
% A=((2*phi1*exp(.5*(phi1+gamma)*T))./((phi1+gamma)*(exp(phi1*T)-1)+2*phi1)).^(2*eta/alpha);
% PriceT=A*exp(-B*r0);
% YieldT=-log(PriceT)./T;
%
% NB: logFlag=1, then
% gamma=exp(vec(1)); eta=exp(vec(2)); alpha=exp(vec(3));

% Check number of inputs ("nargin" counts number of inputs in the function)
	% if there is a LogFlag, and LogFlag=1, then the inputs are in logs
   % Therefore, need to take exp(vec) to obtain the actual parameters
   
if nargin==4
   if logFlag==1
      gamma=exp(vec(1)); eta=exp(vec(2)); alpha=exp(vec(3));
   else
    	gamma=vec(1); eta=vec(2); alpha=vec(3);
	end
elseif nargin<4
	gamma=vec(1); eta=vec(2); alpha=vec(3);
end

% impose restrictions to make sure prices are defined in CIR
% Return NaN if prices are not defined.
if gamma^2+2*alpha<0
   warning('gamma^2+2*alpha<0: Price not Defined')
   PriceT = NaN;
   YieldT = NaN;
elseif eta<= alpha/2
   warning('eta<= alpha/2: Price not Defined')
   PriceT = NaN;
   YieldT = NaN;
else
   
   % if parameters are OK, use the formula to compute
   % the CIR prices.
   
phi1 = sqrt(gamma^2+2*alpha);
B = 2*(exp(phi1*T)-1)./((gamma+phi1)*(exp(phi1*T)-1)+2*phi1);
A =( (2*phi1*exp(.5*(phi1+gamma)*T))./((phi1+gamma)*(exp(phi1*T)-1)+2*phi1)).^(2*eta/alpha);
PriceT = A.*exp(-B*r0);
YieldT = -log(PriceT)./T;

end