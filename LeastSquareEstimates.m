function [L] = LeastSquareEstimates(x,samples,T,TauOne,TauTwo)
Lhood=0; 
for i=1:length(samples)
mu(i)= x(1)*(exp((x(2)-(i)*T)/TauOne) - exp((x(2)-(i)*T)/TauTwo));                       
Lhood = Lhood + (samples(i)-mu(i))^2;
end
L= (Lhood); 
end
