% %% Ordinary least square with penalisation co-efficients.
% %make them zero for resutls with ordinary least square. 
% function [MAP] = MaxAposterioriProb(x,samples,T)
% LstSqrFit=0; lamda=0;lamda_prime=0;
% for i=1:length(samples)
% mu(i)= x(1)*(exp((x(2)-i*T)/x(3))- exp((x(2)-i*T)/x(4))); %no pairwise difference of the samples
% LstSqrFit = LstSqrFit + ((samples(i)-mu(i))^2+lamda*(x(1)^2)+lamda_prime*(5-(x(3)/x(4))^2));
% end
% MAP= (LstSqrFit); 
% end

%% Ordinary least square with priors(MAP)
function [MAP] = MaxAposterioriProb(x,samples,T)
LstSqrFit=0;
for i=1:length(samples)
mu(i)= x(1)*(exp((x(2)-i*T)/x(3))- exp((x(2)-i*T)/x(4))); %no pairwise difference of the samples
LstSqrFit = LstSqrFit + ((samples(i)-mu(i))^2);
end
%% Prior for alpha
a=3; beta=1;
PriorAlpha=(x(1)^(a-1))*exp(-x(1)/beta);
%% Prior for TauOne
a_prime = 1.5; b_prime=1;
PriorTauOne = (x(3)^(a_prime-1))*exp(-x(3)/b_prime);
%% Prior for ratio of TauTwo/TauOne (beta)
a_tilda=3; b_tilda=2;
PriorRatio=((x(4)/x(3))^(a_tilda-1))*((1-(x(4)/x(3)))^(b_tilda-1));
%MAP= -(-(LstSqrFit)+log(PriorAlpha)+log(PriorTauOne)+log(PriorRatio));
MAP= -(exp(-(LstSqrFit))*(PriorAlpha)*(PriorTauOne)*(PriorRatio));
end
