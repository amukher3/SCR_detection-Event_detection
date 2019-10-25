function [MAP] = MaxAposterioriProb(x,samples,T)
% function [MAP] = MaxAposterioriProb(x,samples,T,lambda,lambda_prime,a,beta,a_prime,b_prime,a_tilda,b_tilda)
Lhood=0; 
lambda=1e100;lambda_prime = 1e50;
for i=1:length(samples)
mu(i)= x(1)*(exp((x(2)-i*T)/x(3))- exp((x(2)-i*T)/x(4))); %no pairwise difference of the samples
Lhood = Lhood + ((samples(i)-mu(i))^2+lambda*(x(1)^2)+lambda_prime*(1-(x(3)/x(4))^2));
end
%% Prior for alpha
a=3; beta=1;
PriorAlpha= (x(1)^(a-1))*exp(-x(1)/beta);
%% Prior for TauOne
a_prime =1.5; b_prime=1;
PriorTauOne = (x(3)^(a_prime-1))*exp(-x(3)/b_prime);
%% Prior for ratio of TauTwo/TauOne (beta)
a_tilda=3; b_tilda=2;

PriorRatio=((x(3)/x(4))^(a_tilda-1))*((1-(x(3)/x(4)))^(b_tilda-1));

MAP= (-1*Lhood*PriorAlpha*PriorTauOne*PriorRatio); 
% MAP= (Lhood); 