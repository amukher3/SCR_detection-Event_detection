
function [L,g] = LikelihoodUnderAlternativeWithDerivative(x,sample_differences,noise_var,T)
 Lhood=0; dL_dalpha=0; dL_dtau=0; dL_dtauOne=0; dL_dtauTwo=0;
for i=1:length(sample_differences)
mu(i)= x(1)*(exp((x(2)-(i)*T)/x(3))- exp((x(2)-(i)*T)/x(4))); %no pairwise difference of the samples
Lhood = Lhood + (-1/2*log(2*pi*(noise_var))-(((sample_differences(i)-mu(i))^2)/(2*(noise_var))));

%% derivative part
dmu_dalpha= 1*(exp((x(2)-i*T)/x(3))-exp((x(2)-i*T)/x(4)));   %dmu/d(alpha)
dmu_dtau = x(1)*((exp((x(2)-i*T)/x(3)))/x(3)-(exp((x(2)-i*T)/x(4)))/x(4));    %dmu/d(tau)
dmu_dtauOne = -x(1)*(exp((x(2)-i*T)/x(3))*((x(2)-i*T)/(x(3)^2)));    %dmu/d(tauOne)
dmu_dtauTwo= x(1)*(exp((x(2)-i*T)/x(4))*((x(2)-i*T)/(x(4)^2)));      %dmu/d(tauTwo)
dL_dmu= (sample_differences(i)- mu(i))/(2*(noise_var)^2);            %dL/d(mu)
dL_dalpha = dL_dalpha + dL_dmu*dmu_dalpha;                           %dL/d(alpha)= dL/dmu*dmu/d(alpha)
dL_dtau = dL_dtau + dL_dmu*dmu_dtau;                                 %dL/d(tau)= dL/dmu*dmu/d(tau)
dL_dtauOne = dL_dtauOne + dL_dmu*dmu_dtauOne;                        %dL/d(tauOne)= dL/dmu*dmu/d(tauOne)
dL_dtauTwo = dL_dtauTwo + dL_dmu*dmu_dtauTwo;                        %dL/d(tauTwo)= dL/dmu*dmu/d(tauTwo)

end
L=-(Lhood); 
L_prime_alpha = (dL_dalpha); L_prime_tau =  (dL_dtau); L_prime_tauOne = (dL_dtauOne); L_prime_tauTwo = (dL_dtauTwo);
g=[-L_prime_alpha; -L_prime_tau; -L_prime_tauOne; -L_prime_tauTwo];
end