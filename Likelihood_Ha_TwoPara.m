%function [L,g] = Likelihood_Ha_TwoPara(x,sample_differences,noise_var,T,TauOne,TauTwo)
function [L] = Likelihood_Ha_TwoPara(x,sample_differences,noise_var,T,TauOne,TauTwo)
 Lhood=0; dL_dalpha=0; dL_dtau=0; lambda=0;
for i=1:length(sample_differences)
mu(i)= x(1)*(exp((x(2)-(i)*T)/TauOne)- exp((x(2)-(i)*T)/TauTwo)); %no pairwise difference of the samples
Lhood = Lhood + (((sample_differences(i)-mu(i))^2)+lambda*(x(1))^2);

%% derivative part
% dmu_dalpha= 1*(exp((x(2)-i*T)/TauOne)- exp((x(2)-i*T)/TauTwo));   %dmu/d(alpha)
% dmu_dtau = x(1)*((exp((x(2)-i*T)/TauOne))/TauOne - (exp((x(2)-i*T)/TauTwo))/TauTwo);    %dmu/d(tau)
% dL_dmu = (sample_differences(i)- mu(i))/((noise_var));               %dL/d(mu)
% dL_dalpha = dL_dalpha + dL_dmu*dmu_dalpha;                           %dL/d(alpha)= dL/dmu*dmu/d(alpha)
% dL_dtau = dL_dtau + dL_dmu*dmu_dtau;                                 %dL/d(tau)= dL/dmu*dmu/d(tau)

end
L=(Lhood); 
% L_prime_alpha = (dL_dalpha); L_prime_tau =  (dL_dtau);
% g=[-L_prime_alpha; -L_prime_tau;];
end
