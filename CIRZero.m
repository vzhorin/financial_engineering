function [J,ZZhat]=CIR(vec,Price,Maturity,r0);

global JC_1

gamma=exp(vec(1)); eta=exp(vec(2)); alpha=exp(vec(3)); beta=0;

if gamma^2+2*alpha<=0
   J=JC_1*50;
elseif alpha<=0
   J=JC_1*50;
elseif eta<= - beta*gamma/alpha+alpha/2
   J=JC_1*50;
else
   
Tau=Maturity;

ZZhat=CIRFun([gamma,eta,alpha],r0,Maturity);


if min(isreal(ZZhat))==0
   J=JC_1*50;
else
J=sum(((Price - ZZhat)).^2);

end



end
JC_1=J;